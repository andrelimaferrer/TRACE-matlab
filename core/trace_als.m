function est = trace_als(Y,R,opts)
%TRACE_ALS Step 2 of Algorithm 1: multi-start ALS for the PARAFAC fit.
%
%   est = TRACE_ALS(Y,R)
%   est = TRACE_ALS(Y,R,opts)
%
% Solves the tensor-fitting problem of eq. (17),
%
%   min_{T_A,T_H,T_S} || Y - I_{3,L} x_1 T_A x_2 T_H x_3 T_S ||_F^2,
%
% by alternating least squares. Each of the three mode unfoldings
% (eqs. (18)-(20)) is linear in one factor, so each factor admits the
% closed-form conditional update of eq. (21),
%
%   T_X_hat = Y_(n) (Z_X_hat^T)^dagger,      X in {A,H,S}, n = 1,2,3,
%
%   Z_A = T_S <> T_H,   Z_H = T_S <> T_A,   Z_S = T_H <> T_A.
%
% Each update is exact given the other two factors, and the cycle needs
% no external optimization package.
%
% THREE IMPLEMENTATION CHOICES THAT MATTER
%
%   1. SVD INITIALIZATION (Algorithm 1, Step 1). The first start takes the
%      L dominant left singular vectors of the corresponding unfolding --
%      a consistent estimate of each factor's column space in the
%      noiseless case -- instead of a random guess.
%
%   2. QR/BACKSLASH UPDATES. Eq. (21) is written with an explicit
%      pseudo-inverse, Y_(n) Z^* (Z^T Z^*)^{-1}. Forming the Gram matrix
%      Z^T Z^* SQUARES the condition number, which is ruinous precisely in
%      the ill-conditioned regime the paper studies (the T = 4 point of
%      Fig. (c) has EM coherence 0.985). The code therefore solves the LS
%      problem directly with MATLAB's backslash, which uses a QR
%      factorization of Z. The result is the same least-squares solution,
%      computed stably.
%
%   3. MULTI-START. Problem (17) is non-convex. The cycle is repeated from
%      N_s - 1 further INDEPENDENT random initializations and the
%      smallest-residual fit is retained. The paper uses N_s = 5,
%      tolerance eps = 1e-8 and I_max = 200.
%
% Between sweeps the factors are renormalized (TRACE_NORMALIZE_FACTORS)
% to stop the CP scale from drifting; this changes no represented tensor.
%
% ON THE CHOICE OF SOLVER (Remark 2 of the paper)
%   ALS is a block-coordinate method and converges only linearly.
%   Gauss-Newton solvers treat the three factors jointly and are the
%   standard remedy for the swamps and bottlenecks that arise when two
%   components become nearly collinear in one mode. TRACE is agnostic:
%   both solve (17) and are followed by the same recovery stage. ALS is
%   adopted here because it is exact per block and, at the operating
%   points of Section V, was found indistinguishable from a Gauss-Newton
%   solver. This package ships the ALS path only, which is why it has no
%   third-party dependency.
%
% MODEL CONVENTION
%   Y in C^(NRF x T x P) = [[T_A,T_H,T_S]], eq. (11).
%
% INPUTS
%   Y                 received tensor, NRF x T x P
%   R                 CP rank (R = L for single-user TRACE, eq. (15))
%
%   opts.maxit        maximum ALS sweeps I_max (default 200)
%   opts.tol          stopping tolerance epsilon (default 1e-8), applied
%                     to the RELATIVE DECREASE of the residual
%   opts.verbose      print progress (default false)
%   opts.nstarts      number of starts N_s (default 1; paper uses 5)
%   opts.init_method  first-start initialization:
%                       'svd'      (default, Algorithm 1 Step 1)
%                       'random'
%                       'provided'
%   opts.init         supplied factors for init_method='provided'
%
% OUTPUTS
%   est.TA, est.TH, est.TS   factor estimates, eqs. (12)-(14)
%   est.fit           residual history for the selected start
%   est.tensor_nmse   final tensor fitting NMSE (the objective of (17),
%                     normalized by ||Y||^2) -- this is the tensor fit,
%                     NOT the channel NMSE reported in Section V
%   est.runtime       total runtime over all starts
%   est.best_runtime  runtime of the selected start
%   est.iterations    iterations of the selected start
%   est.nstarts       number of starts used
%   est.best_start    selected start index
%   est.start_nmse    final tensor NMSE of every start
%   est.method        'TRACE-ALS-SVD-QR-Multistart'
%
% COMPLEXITY
%   O(P T NRF L + L^3) per sweep, N_s times that per estimate
%   (Section VI, complexity paragraph).
%
% PAPER CROSS-REFERENCE
%   eqs. (17)-(21); Algorithm 1, Steps 1-2; Remark 2
%
% See also TRACE_PHYSICAL_RECOVERY, TRACE_UNFOLD, TRACE_NORMALIZE_FACTORS.

if nargin < 3 || isempty(opts)
    opts = struct();
end

if ~isfield(opts,'maxit'),       opts.maxit = 200; end
if ~isfield(opts,'tol'),         opts.tol = 1e-8; end
if ~isfield(opts,'verbose'),     opts.verbose = false; end
if ~isfield(opts,'nstarts'),     opts.nstarts = 1; end
if ~isfield(opts,'init_method'), opts.init_method = 'svd'; end
if ~isfield(opts,'init'),        opts.init = []; end

if opts.nstarts < 1 || opts.nstarts ~= floor(opts.nstarts)
    error('opts.nstarts must be a positive integer.');
end

%% Precompute the unfoldings (18)-(20) once and reuse across all starts.
Y1 = trace_unfold(Y,1);   % NRF x (T P)   eq. (18)
Y2 = trace_unfold(Y,2);   % T   x (NRF P) eq. (19)
Y3 = trace_unfold(Y,3);   % P   x (NRF T) eq. (20)

bestNMSE = Inf;
bestEst = [];
bestStart = 0;
start_nmse = Inf(opts.nstarts,1);
totalRuntime = 0;

for s = 1:opts.nstarts

    localopts = opts;

    % The requested initialization is used for the first start. All
    % additional starts are genuinely independent random initializations.
    if s == 1
        localopts.init_method = opts.init_method;
        localopts.init = opts.init;
    else
        localopts.init_method = 'random';
        localopts.init = [];
    end

    est_s = trace_als_single(Y,Y1,Y2,Y3,R,localopts);

    start_nmse(s) = est_s.tensor_nmse;
    totalRuntime = totalRuntime + est_s.runtime;

    if opts.verbose
        fprintf('ALS start %d/%d: tensor NMSE = %.4e\n', ...
            s,opts.nstarts,est_s.tensor_nmse);
    end

    % Keep the smallest residual, as prescribed by Algorithm 1 Step 2.
    if est_s.tensor_nmse < bestNMSE
        bestNMSE = est_s.tensor_nmse;
        bestEst = est_s;
        bestStart = s;
    end
end

est = bestEst;
est.runtime = totalRuntime;
est.best_runtime = bestEst.runtime;
est.nstarts = opts.nstarts;
est.best_start = bestStart;
est.start_nmse = start_nmse;
est.method = 'TRACE-ALS-SVD-QR-Multistart';

end


%% =====================================================================
function est = trace_als_single(Y,Y1,Y2,Y3,R,opts)
%TRACE_ALS_SINGLE One ALS run from one initialization.

%% Initialization
switch lower(opts.init_method)

    case 'svd'
        [TA,TH,TS] = svd_initialization(Y1,Y2,Y3,Y,R);

    case 'random'
        init = trace_random_init(Y,R);
        TA = init.TA;
        TH = init.TH;
        TS = init.TS;

    case 'provided'
        if isempty(opts.init)
            error(['opts.init_method is ''provided'', but opts.init ', ...
                   'was not supplied.']);
        end
        TA = opts.init.TA;
        TH = opts.init.TH;
        TS = opts.init.TS;

    otherwise
        error(['Unknown opts.init_method. Use ''svd'', ''random'', ', ...
               'or ''provided''.']);
end

[TA,TH,TS] = trace_normalize_factors(TA,TH,TS);

fit = zeros(opts.maxit+1,1);
fit(1) = relative_error(Y,TA,TH,TS);

if opts.verbose
    fprintf('  initialization (%s): rel. error = %.4e\n', ...
        opts.init_method,fit(1));
end

%% ALS iterations: one sweep = the three conditional updates of eq. (21)
tic;
nit = opts.maxit;

for it = 1:opts.maxit

    % Mode-1 update, eq. (18): Y_(1) = T_A (T_S <> T_H)^T.
    % (ZA \ Y1.').' is the LS solution of min ||Y1 - TA*ZA.'||_F, solved
    % by QR rather than through the normal equations of eq. (21).
    ZA = khatrirao(TS,TH);                 % (P*T) x R
    TA = (ZA \ Y1.').';                    % direct LS via MATLAB backslash

    % Mode-2 update, eq. (19): Y_(2) = T_H (T_S <> T_A)^T.
    ZH = khatrirao(TS,TA);                 % (P*NRF) x R
    TH = (ZH \ Y2.').';

    % Mode-3 update, eq. (20): Y_(3) = T_S (T_H <> T_A)^T.
    ZS = khatrirao(TH,TA);                 % (T*NRF) x R
    TS = (ZS \ Y3.').';

    % Control CP scale drift while preserving the represented tensor.
    [TA,TH,TS] = trace_normalize_factors(TA,TH,TS);

    fit(it+1) = relative_error(Y,TA,TH,TS);

    if opts.verbose && (it == 1 || mod(it,10) == 0)
        fprintf('  ALS iter %3d: rel. error = %.4e\n',it,fit(it+1));
    end

    % Stop on the relative DECREASE, not on the residual level: the noise
    % floor sets the attainable residual, so an absolute threshold would
    % be SNR-dependent.
    relDecrease = abs(fit(it)-fit(it+1))/max(fit(it),eps);

    if relDecrease <= opts.tol
        nit = it;
        fit = fit(1:it+1);
        break;
    end
end

elapsed = toc;

Yhat = trace_cpdgen(TA,TH,TS);
tensor_nmse = norm(Y(:)-Yhat(:))^2/max(norm(Y(:))^2,eps);

est = struct();
est.TA = TA;
est.TH = TH;
est.TS = TS;
est.fit = fit;
est.tensor_nmse = tensor_nmse;
est.runtime = elapsed;
est.iterations = nit;
est.init_method = opts.init_method;

end


%% =====================================================================
function [TA,TH,TS] = svd_initialization(Y1,Y2,Y3,Y,R)
%SVD_INITIALIZATION Algorithm 1, Step 1.
%
% Takes the R dominant left singular vectors of the mode-1 and mode-2
% unfoldings as initial T_A and T_H. In the noiseless case the column
% space of Y_(n) equals the column space of the corresponding factor
% (eqs. (18)-(19)), so this is a consistent estimate of each factor's
% subspace -- it does not identify the individual columns, which is what
% the ALS cycle then does.
%
% T_S is not taken from a third SVD: with T_A and T_H in hand, eq. (20)
% gives it exactly by one least-squares solve, which is both cheaper and
% consistent with the other two.

NRF = size(Y1,1);
T = size(Y2,1);

% RF-mode dominant subspace.
[UA,~,~] = svd(Y1,'econ');
qA = min(R,size(UA,2));
TA = zeros(NRF,R,'like',Y);
TA(:,1:qA) = UA(:,1:qA);
% Degenerate case NRF < R: pad with random columns so the factor is full
% width. This cannot happen when eq. (29) is satisfied.
if qA < R
    TA(:,qA+1:R) = randn(NRF,R-qA)+1j*randn(NRF,R-qA);
end

% EM-mode dominant subspace.
[UH,~,~] = svd(Y2,'econ');
qH = min(R,size(UH,2));
TH = zeros(T,R,'like',Y);
TH(:,1:qH) = UH(:,1:qH);
if qH < R
    TH(:,qH+1:R) = randn(T,R-qH)+1j*randn(T,R-qH);
end

% Estimate T_S by direct LS from the mode-3 unfolding, eq. (20).
ZS = khatrirao(TH,TA);
TS = (ZS \ Y3.').';

% Match the initialization energy to the data.
Y0 = trace_cpdgen(TA,TH,TS);
if norm(Y0(:)) > eps
    TS = TS*(norm(Y(:))/norm(Y0(:)));
end

end


%% =====================================================================
function e = relative_error(Y,TA,TH,TS)
%RELATIVE_ERROR ||Y - [[TA,TH,TS]||_F / ||Y||_F, the objective of eq. (17).
Yhat = trace_cpdgen(TA,TH,TS);
e = norm(Y(:)-Yhat(:))/max(norm(Y(:)),eps);
end

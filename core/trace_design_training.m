function [training,diag] = trace_design_training(par)
%TRACE_DESIGN_TRAINING Design the EM probing and RF combining matrices.
%
%   [training,diag] = TRACE_DESIGN_TRAINING(par)
%
% Builds the two known training matrices of the paper,
%
%   F_EMbar = [fbar_{EM,1},...,fbar_{EM,T}]^T  in C^(T x M)   (eq. (9))
%   F_RF                                        in C^(NRF x Nr)
%
% and returns a diagnostics structure that certifies the design is not
% degenerate. Both matrices are generated ONCE per experiment and reused
% for every Monte-Carlo realization and every SNR point, which reflects
% the fixed training protocol assumed in Section V.
%
% WHY RANDOM PHASES ARE THE DEFAULT
%   F_EMbar characterizes the array: designing it means choosing which M
%   reconfiguration states to cycle through. The unit-modulus,
%   independent-uniform-phase model is a statistical surrogate for a
%   sufficiently diverse set of states (see the discussion following
%   eq. (16)), adopted because it yields low mutual coherence
%
%       mu(F_EMbar) = max_{m != m'} |fbar_m^H fbar_m'|
%                                   / (||fbar_m|| ||fbar_m'||)   (eq. (16))
%
%   on which the support detection of eq. (23) relies. F_RF, by contrast,
%   is a commanded weight matrix realized by phase shifters, so
%   constant-modulus entries are a hardware constraint, not a modelling
%   choice.
%
% THE DFT TRAP (Section VI of the paper)
%   The 'partial-dft' option is provided for completeness but is NOT the
%   design used for the published curves. Subsampling an M-point DFT at
%   the EQUALLY SPACED rows floor((t-1)M/T)+1 makes columns m and m+T of
%   F_EMbar identical whenever M/T is an integer -- the omitted rows are
%   exactly the ones that would distinguish them. The Kruskal rank k_H
%   then collapses to 1 and the sufficient condition (28) fails. For
%   T = 12, M = 24 this produces 12 duplicate pairs. To avoid handing the
%   reader a booby-trapped default, the partial-DFT branch below selects
%   DISTINCT RANDOM rows (randperm), never equally spaced ones. The
%   duplicate-column detector further catches the failure if it occurs.
%   Deterministic codebooks are safe when T >= M.
%
% INPUT fields used
%   par.T, par.M, par.NRF, par.Nr
%   par.FEMbar_type : 'random-phase' (default) or 'partial-dft'
%   par.FRF_type    : 'random-phase' (default) or 'partial-dft'
%
% OUTPUTS
%   training.FEMbar : T x M EM probing matrix, eq. (9)
%   training.FRF    : NRF x Nr RF combining matrix
%   training.diag   : diagnostics structure (same as second output)
%
%   diag.mu_EM                  mutual coherence of eq. (16)
%   diag.duplicate_EM_pairs     number of essentially duplicate EM columns
%   diag.rank_FEMbar            numerical rank of F_EMbar
%   diag.rank_FRF               numerical rank of F_RF
%   diag.min_RF_gain_rel        minimum RF angular gain / maximum gain
%   diag.theta_min_RF_gain      angle of the minimum RF gain (radians)
%   diag.is_degenerate          true if a severe issue is detected
%
% THE TWO DIAGNOSTICS AND WHAT THEY PROTECT
%   duplicate_EM_pairs > 0 means k_H = 1 and the uniqueness condition (28)
%   is void; min_RF_gain_rel near zero means F_RF has a deep angular null,
%   so paths arriving near theta_min_RF_gain are invisible to the array
%   and the corresponding column of T_A = F_RF A(theta) (eq. (12)) is
%   nearly zero. Both conditions make the recovered channel meaningless
%   while leaving the tensor fit apparently healthy, so they are reported
%   loudly rather than silently tolerated.
%
% PAPER CROSS-REFERENCE
%   eqs. (9), (16), (28), (29); Section VI, degenerate-codebook discussion
%
% See also TRACE_GENERATE_CLEAN_DATA, TRACE_PREPARE_RECOVERY.

T   = par.T;
M   = par.M;
NRF = par.NRF;
Nr  = par.Nr;

if ~isfield(par,'FEMbar_type'), par.FEMbar_type = 'random-phase'; end
if ~isfield(par,'FRF_type'),    par.FRF_type    = 'random-phase'; end

%% ---------------------------------------------------------
% EM probing matrix F_EMbar, eq. (9): T x M, one row per probing state
%% ---------------------------------------------------------
switch lower(par.FEMbar_type)

    case 'random-phase'
        % Constant-modulus complex probing. The 1/sqrt(T) scaling gives
        % every column unit norm, so mu(F_EMbar) in eq. (16) is read
        % directly off the Gram matrix and the SNR definition of
        % Section V is independent of T.
        FEMbar = exp(1j*2*pi*rand(T,M))/sqrt(T);

    case 'partial-dft'
        if T >= M
            % Tall DFT frame with orthonormal columns: mu = 0, safe.
            D = local_dftmat(T)/sqrt(T);
            FEMbar = D(:,1:M);
        else
            % Select DISTINCT nonperiodic rows. Random row selection avoids
            % the deterministic aliasing caused by equally spaced rows
            % (see THE DFT TRAP above).
            D = local_dftmat(M)/sqrt(M);
            rows = randperm(M,T);
            FEMbar = D(rows,:);
        end

    otherwise
        error('Unknown par.FEMbar_type. Use random-phase or partial-dft.');
end

%% ---------------------------------------------------------
% RF combining matrix F_RF: NRF x Nr, realized by phase shifters
%% ---------------------------------------------------------
switch lower(par.FRF_type)

    case 'random-phase'
        % Constant-modulus phase-shifter combiner. Each row has unit norm.
        FRF = exp(1j*2*pi*rand(NRF,Nr))/sqrt(Nr);

    case 'partial-dft'
        if NRF >= Nr
            D = local_dftmat(NRF)/sqrt(NRF);
            FRF = D(:,1:Nr);
        else
            % Select DISTINCT DFT rows rather than periodic/equispaced rows.
            D = local_dftmat(Nr)/sqrt(Nr);
            rows = randperm(Nr,NRF);
            FRF = D(rows,:);
        end

    otherwise
        error('Unknown par.FRF_type. Use random-phase or partial-dft.');
end

%% ---------------------------------------------------------
% Training-design diagnostics
%% ---------------------------------------------------------
% Mutual coherence of eq. (16): normalized EM column Gram matrix with the
% diagonal blanked out, then the largest off-diagonal magnitude.
Fn = FEMbar ./ max(vecnorm(FEMbar),eps);
GEM = abs(Fn'*Fn);
GEM(1:size(GEM,1)+1:end) = 0;

if isempty(GEM)
    mu_EM = 0;
else
    mu_EM = max(GEM(:));
end

% Coherence exactly 1 means two EM columns are collinear: the two
% corresponding reconfiguration states are indistinguishable and k_H
% collapses, voiding eq. (28).
duplicate_threshold = 1 - 1e-10;
duplicate_EM_pairs = nnz(triu(GEM > duplicate_threshold,1));

% Numerical ranks.
rank_FEMbar = rank(FEMbar);
rank_FRF = rank(FRF);

% Scan the RF-compressed steering gain ||F_RF a(theta)||^2 over the
% physical angular range. A deep null here means some AoAs produce a
% near-zero column of T_A (eq. (12)) and cannot be estimated at all.
theta_diag = linspace(-pi/2,pi/2,4001);
A_diag = ula_steering(Nr,theta_diag);
rf_gain = sum(abs(FRF*A_diag).^2,1);

[max_RF_gain,~] = max(rf_gain);
[min_RF_gain,idx_min] = min(rf_gain);
min_RF_gain_rel = min_RF_gain/max(max_RF_gain,eps);
theta_min_RF_gain = theta_diag(idx_min);

% Severe-degeneracy tests. These thresholds are intentionally conservative:
% they flag exact/near-exact duplicate EM columns or a very deep RF null.
deep_null_threshold = 1e-3;
is_degenerate = (duplicate_EM_pairs > 0) || ...
                (min_RF_gain_rel < deep_null_threshold);

diag = struct();
diag.mu_EM = mu_EM;
diag.duplicate_EM_pairs = duplicate_EM_pairs;
diag.rank_FEMbar = rank_FEMbar;
diag.rank_FRF = rank_FRF;
diag.min_RF_gain_rel = min_RF_gain_rel;
diag.theta_min_RF_gain = theta_min_RF_gain;
diag.is_degenerate = is_degenerate;
diag.FEMbar_type = par.FEMbar_type;
diag.FRF_type = par.FRF_type;

%% Warnings: never allow a degenerate codebook to pass silently.
if duplicate_EM_pairs > 0
    warning('TRACE:DegenerateEMCodebook', ...
        ['Degenerate EM codebook: %d essentially duplicate column pair(s) ', ...
         'detected (max coherence %.6g). Change the training design.'], ...
        duplicate_EM_pairs,mu_EM);
end

if min_RF_gain_rel < deep_null_threshold
    warning('TRACE:DeepRFNull', ...
        ['RF training matrix has a deep angular null: minimum normalized ', ...
         'gain %.3e at theta = %.2f deg. Change the training design.'], ...
        min_RF_gain_rel,theta_min_RF_gain*180/pi);
end

training = struct();
training.FEMbar = FEMbar;
training.FRF = FRF;
training.diag = diag;

end

%% =====================================================================
function D = local_dftmat(N)
%LOCAL_DFTMAT Unnormalized N-point DFT matrix.
% Equivalent to dftmtx(N) from the Signal Processing Toolbox, written out
% so that this package needs no toolbox at all.
n = (0:N-1);
D = exp(-2j*pi*(n.')*n/N);
end

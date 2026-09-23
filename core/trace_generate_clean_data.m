function data = trace_generate_clean_data(par,training)
%TRACE_GENERATE_CLEAN_DATA One noiseless realization of the TRACE model.
%
%   data = TRACE_GENERATE_CLEAN_DATA(par,training)
%
% Draws one channel realization (AoAs, EM signatures, path gains) and
% synthesizes the noiseless pilot tensor of eq. (11). For the single-user
% uplink of the paper K = 1 and the tensor rank equals the number of
% resolvable paths, R = L (eq. (15)).
%
% MODEL IMPLEMENTED
%   Y = [[T_A, T_H, T_S]]  in C^(NRF x T x P)                   eq. (11)
%     T_A = F_RF A(theta)                                       eq. (12)
%     T_H = F_EMbar H_EMbar                                     eq. (13)
%     T_S = B = diag(s) C                                       eqs. (10),(14)
%
%   with s in C^P the known pilot sequence and C in C^(P x L) the stacked
%   path-gain matrix (called `alpha` in this code, see the symbol table in
%   README.md).
%
% RANDOM DRAWS, MATCHING SECTION V
%   AoAs      : uniform in [theta_min,theta_max] = [-60,60] degrees,
%               rejection-sampled to enforce a minimum separation
%               (15 degrees for Fig. (a), 5 degrees for Figs. (b) and (c))
%   EM states : one per path, drawn WITHOUT replacement from the M states,
%               giving the one-hot columns of H_EMbar in eq. (3)
%   gains     : i.i.d. CN(0,1/L) across paths and slots. Note the 1/L:
%               the code line (randn+1j*randn)/sqrt(2*L) has variance
%               1/(2L) per real part, hence 1/L total, so the aggregate
%               channel power is independent of L
%   pilot     : s_p = 1/sqrt(P), constant modulus with |s_p|^2 = 1/P and
%               unit total energy, as stated below eq. (4)
%
% WHY THE ANGLE SEPARATION IS ENFORCED
%   Two AoAs closer than the array's Rayleigh resolution make two columns
%   of T_A nearly collinear. Kruskal's condition (28) is then satisfied
%   only nominally, the CP problem becomes a bottleneck, and the reported
%   NMSE would measure the pathology rather than the estimator. The
%   minimum separation is part of the experimental setup, not a fix.
%
% INPUT fields
%   par.L, par.Nr, par.NRF, par.M, par.T, par.P
%   optional: par.theta_min (-pi/3), par.theta_max (pi/3),
%             par.min_angle_sep_deg (15)
%   training : output of TRACE_DESIGN_TRAINING (F_EMbar and F_RF)
%
% OUTPUT
%   data.Yclean  noiseless tensor                            eq. (11)
%   data.TA/.TH/.TS   true latent factors                    eqs. (12)-(14)
%   data.A       Nr x L steering matrix A(theta)             eq. (2)
%   data.HEM     M x L one-hot EM signature matrix H_EMbar   eq. (3)
%   data.theta   1 x L true AoAs (radians)
%   data.alpha   P x L true path gains, the matrix C of eq. (10)
%   data.s       P x 1 known pilot sequence
%   data.FEMbar, data.FRF, data.em_support, data.R, data.par
%
% PAPER CROSS-REFERENCE
%   eqs. (2), (3), (9)-(15); Section V, simulation setup paragraph
%
% See also TRACE_DESIGN_TRAINING, TRACE_ADD_NOISE, TRACE_CPDGEN.

L=par.L;
R=L;                     % single user: tensor rank = number of paths
Nr=par.Nr;
M=par.M;
P=par.P;
FRF=training.FRF;
FEMbar=training.FEMbar;

%% One-hot EM angular signatures: H_EMbar of eq. (3), then T_H of eq. (13)
% Each path occupies exactly one of the M EM angular samples. Drawing
% without replacement (R <= M) guarantees distinct EM signatures, which is
% the premise behind k_H = min(T,L) in eq. (29).
HEM=zeros(M,R);
if R<=M
    em_support=randperm(M,R);
else
    em_support=randi(M,1,R);
end
HEM(sub2ind([M R],em_support,1:R))=1;
TH=FEMbar*HEM;

%% Spatial angles and steering matrix: A(theta) of eq. (2), T_A of eq. (12)
if ~isfield(par,'theta_min'), par.theta_min=-pi/3; end
if ~isfield(par,'theta_max'), par.theta_max= pi/3; end
if ~isfield(par,'min_angle_sep_deg'), par.min_angle_sep_deg=15; end

sep=par.min_angle_sep_deg*pi/180;
range=par.theta_max-par.theta_min;

if R>1 && (R-1)*sep>range
    error('Requested minimum angular separation is infeasible for L paths.');
end

% Rejection sampling: draw uniformly, keep only candidates far enough from
% all previously accepted angles. The attempt cap turns an infeasible
% request into a clear error instead of a hang.
theta=zeros(1,R);
count=0;
attempts=0;

while count<R
    attempts=attempts+1;
    if attempts>1e5
        error('Could not generate L separated angles.');
    end

    cand=par.theta_min+range*rand;
    if count==0 || all(abs(cand-theta(1:count))>=sep)
        count=count+1;
        theta(count)=cand;
    end
end

A=ula_steering(Nr,theta);
TA=FRF*A;

%% Known pilot sequence s, |s_p|^2 = 1/P
% Unit-energy deterministic pilot. For K=1 only a scalar pilot is needed
% at each slot; there is no user-separation pilot matrix.
s=ones(P,1)/sqrt(P);

%% Path gains C and the third factor T_S = B = diag(s) C, eqs. (10),(14)
% CN(0,1/L) per entry: variance 1/(2L) in each real dimension.
alpha=(randn(P,R)+1j*randn(P,R))/sqrt(2*L);
TS=alpha.*s;                    % B = diag(s)*alpha, eq. (10)

%% Noiseless received tensor, eq. (11)
Yclean=trace_cpdgen(TA,TH,TS);

%% Output structure
data=struct();
data.Yclean=Yclean;
data.TA=TA;
data.TH=TH;
data.TS=TS;
data.A=A;
data.HEM=HEM;
data.theta=theta;
data.alpha=alpha;
data.s=s;
data.FEMbar=FEMbar;
data.FRF=FRF;
data.em_support=em_support;
data.R=R;
data.par=par;
end

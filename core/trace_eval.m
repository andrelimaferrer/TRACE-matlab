function metrics = trace_eval(data,est,phys)
%TRACE_EVAL Simulation-only scoring of one TRACE run.
%
%   metrics = TRACE_EVAL(data,est,phys)
%
% Collects every quantity used to build the figures and the diagnostics
% quoted in Section V. Everything here requires ground truth, so this
% function belongs to the experiment harness, not to the estimator.
%
% WHAT IS A RESULT AND WHAT IS A DIAGNOSTIC
%   metrics.channel_nmse is THE reported metric (Section V). It is
%   end-to-end and permutation-invariant, so it needs no path matching.
%
%   Everything else -- the per-factor NMSEs, the per-parameter NMSEs and
%   the EM support error rate -- is per-path and therefore only defined
%   after the CP column permutation has been resolved against the truth.
%   These are reported nowhere in the paper; they exist to diagnose which
%   stage of Algorithm 1 fails when a realization goes wrong.
%
% RESOLVING THE PERMUTATION AND THE SCALING
%   The matching score sums the moduli of the normalized column inner
%   products across all three modes, so a path is matched only if it
%   agrees in the RF, EM and fading modes at once. The best assignment is
%   then found by TRACE_BEST_PERMUTATION. Scale is removed separately by
%   scaled_nmse below, which projects each estimated column onto its true
%   counterpart before measuring the error -- the least-squares optimal
%   complex scalar, which is the right way to quotient out the CP scaling
%   ambiguity.
%
%   This built-in alignment replaces Tensorlab's cpderr; the package has
%   no third-party dependency.
%
% INPUTS
%   data   ground-truth realization (TRACE_GENERATE_CLEAN_DATA + noise)
%   est    factor estimates from TRACE_ALS
%   phys   physical estimates from TRACE_PHYSICAL_RECOVERY
%
% OUTPUT FIELDS
%   channel_nmse             the reported metric, Section V
%   nmse_TA/_TH/_TS          latent factor NMSEs, eqs. (12)-(14)
%   relerr_TA/_TH/_TS        square roots of the above
%   nmse_HEM/_A/_B/_alpha    physical parameter NMSEs, eqs. (3),(2),(10)
%   HEM_support_error_rate   fraction of paths whose EM state (eq. (23))
%                            was mis-detected
%   perm                     resolved CP column permutation
%   rank_ZS/_deficient_ZS/cond_ZS
%                            conditioning of the LS system of eq. (26)
%
% PAPER CROSS-REFERENCE
%   Section V (channel NMSE); eqs. (2), (3), (10), (12)-(14), (23), (26)
%
% See also TRACE_CHANNEL_NMSE, TRACE_BEST_PERMUTATION, TRACE_ALS.

%% Align latent factors (simulation-only per-factor NMSEs)
TA0=normcols(data.TA); TAe=normcols(est.TA);
TH0=normcols(data.TH); THe=normcols(est.TH);
TS0=normcols(data.TS); TSe=normcols(est.TS);

% Agreement must hold in all three modes simultaneously.
Score=abs(TA0'*TAe)+abs(TH0'*THe)+abs(TS0'*TSe);
perm=trace_best_permutation(Score);

metrics.nmse_TA=scaled_nmse(data.TA,est.TA(:,perm));
metrics.nmse_TH=scaled_nmse(data.TH,est.TH(:,perm));
metrics.nmse_TS=scaled_nmse(data.TS,est.TS(:,perm));

metrics.relerr_TA=sqrt(metrics.nmse_TA);
metrics.relerr_TH=sqrt(metrics.nmse_TH);
metrics.relerr_TS=sqrt(metrics.nmse_TS);

metrics.perm=perm;

%% Physical parameter metrics (per path, hence permuted)
HEM=phys.HEM_hat(:,perm);
A=phys.A_hat(:,perm);
B=phys.B_hat(:,perm);
alpha=phys.alpha_hat(:,perm);

metrics.nmse_HEM=norm(data.HEM(:)-HEM(:))^2 / ...
                 max(norm(data.HEM(:))^2,eps);
metrics.nmse_A=norm(data.A(:)-A(:))^2 / ...
               max(norm(data.A(:))^2,eps);
metrics.nmse_B=norm(data.TS(:)-B(:))^2 / ...
               max(norm(data.TS(:))^2,eps);
metrics.nmse_alpha=norm(data.alpha(:)-alpha(:))^2 / ...
                   max(norm(data.alpha(:))^2,eps);

% Both H_EMbar and its estimate are one-hot (eq. (3), eq. (23)), so the
% support is the row index of the single nonzero entry of each column.
[~,true_supp]=max(abs(data.HEM),[],1);
[~,hat_supp]=max(abs(HEM),[],1);
metrics.HEM_support_error_rate=mean(true_supp~=hat_supp);

%% End-to-end channel NMSE: no oracle permutation is needed
metrics.channel_nmse=trace_channel_nmse( ...
    data.A,data.HEM,data.alpha, ...
    phys.A_hat,phys.HEM_hat,phys.alpha_hat);

%% Recovery conditioning diagnostics (eq. (26))
metrics.rank_ZS=phys.rank_ZS;
metrics.rank_deficient_ZS=phys.rank_deficient_ZS;
metrics.cond_ZS=phys.cond_ZS;
end

%% ---------------------------------------------------------------------
function X=normcols(X)
%NORMCOLS Scale every column to unit norm.
X=X./max(vecnorm(X),eps);
end

%% ---------------------------------------------------------------------
function e=scaled_nmse(A,B)
%SCALED_NMSE Column-wise NMSE after removing the CP scaling ambiguity.
% beta is the least-squares optimal complex scalar per column, so the
% error measured is the part that no rescaling of B could remove.
beta=sum(conj(B).*A,1)./max(sum(abs(B).^2,1),eps);
B=B.*beta;
e=norm(A(:)-B(:))^2/max(norm(A(:))^2,eps);
end

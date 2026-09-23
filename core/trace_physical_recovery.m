function phys = trace_physical_recovery(est,data,cache)
%TRACE_PHYSICAL_RECOVERY Steps 3-7 of Algorithm 1: factors -> physical parameters.
%
%   phys = TRACE_PHYSICAL_RECOVERY(est,data,cache)
%
% Second stage of TRACE. The tensor-fitting stage (TRACE_ALS, Step 2 of
% Algorithm 1) returns factors that satisfy, up to noise, eq. (22):
%
%     T_H_hat ~ F_EMbar H_EMbar,      T_A_hat ~ F_RF A(theta),
%
% so each physical quantity is recovered by matching against its known
% training dictionary. This function implements, in order:
%
%   Step 3  EM support detection            P_EM, eq. (23)
%   Step 4  AoA beamspace search            P_A,  eq. (24)
%   Step 5  structured factors T_H~, T_A~   eq. (25)
%   Step 6  re-estimation of B by LS        eq. (26)
%   Step 7  path gains C = diag(s)^-1 B     eq. (27)
%
% WHY STEP 6 EXISTS
%   Re-estimating B from the STRUCTURED factors, rather than keeping the
%   T_S delivered by ALS, re-references the third factor to the true
%   training matrices and thereby removes the CP scaling ambiguity. This
%   is not a refinement of convenience: without it the recovered gains
%   carry an arbitrary per-column scalar and the channel NMSE is
%   undefined.
%
% PERMUTATION
%   The CP column permutation simply relabels propagation paths. Steps 3-7
%   act column-wise and are permutation-equivariant, so no labelling is
%   required (Section VI, first paragraph). Because K = 1, no user
%   association and no block-diagonal reconstruction are needed either.
%
% INPUTS
%   est    factor estimates from TRACE_ALS (.TA, .TH, .TS)
%   data   realization structure; supplies F_EMbar, F_RF, Y, s, R, par
%   cache  dictionaries from TRACE_PREPARE_RECOVERY
%
% OUTPUTS
%   phys.HEM_hat    M x L, one-hot, eq. (23)
%   phys.A_hat      Nr x L, eq. (24)
%   phys.theta_hat  1 x L estimated AoAs (radians)
%   phys.TH_phys, phys.TA_phys   structured factors of eq. (25)
%   phys.B_hat      P x L, eq. (26)
%   phys.alpha_hat  P x L path gains, the matrix C of eq. (27)
%   phys.rank_ZS, .rank_deficient_ZS, .cond_ZS, .singular_values_ZS
%                   conditioning diagnostics for the LS system of eq. (26)
%
% PAPER CROSS-REFERENCE
%   eqs. (22)-(27); Algorithm 1, Steps 3-7
%
% See also TRACE_ALS, TRACE_PREPARE_RECOVERY, TRACE_EVAL.

R=data.R;
P=data.par.P;

%% Step 3: EM support detection, eq. (23)
% mhat_r = argmax_m |fbar_{EM,m}^H t_{H,r}|, independently for each path.
% The result is one-hot by construction, matching the structure of
% H_EMbar in eq. (3).
scoreH=abs(cache.FEMnorm'*est.TH);
[~,mhat]=max(scoreH,[],1);

HEM_hat=zeros(size(data.HEM),'like',est.TH);
HEM_hat(sub2ind(size(HEM_hat),mhat,1:R))=1;

%% Step 4: AoA recovery in compressed RF beamspace, eq. (24)
% theta_r = argmax_{theta in Theta} |(F_RF a(theta))^H t_{A,r}|.
% The search is over the RF-PROJECTED dictionary, not over a(theta)
% itself: only NRF < Nr projections are observed, so the matched filter
% must live in beamspace. Resolution is limited by the grid Theta, which
% is the source of the error floor seen at high SNR in Fig. (b).
scoreA=abs(cache.PhiRFnorm'*est.TA);
[~,gidx]=max(scoreA,[],1);

theta_hat=cache.theta_grid(gidx);
A_hat=cache.Agrid(:,gidx);

%% Step 5: rebuild the physically scaled factors, eq. (25)
TH_phys=data.FEMbar*HEM_hat;
TA_phys=data.FRF*A_hat;

%% Step 6: re-estimate B = T_S from the mode-3 unfolding, eq. (26)
% B_hat = Y_(3) [(T_H~ <> T_A~)^T]^dagger.
Y3=trace_unfold(data.Y,3);
ZS=khatrirao(TH_phys,TA_phys);             % (T*NRF) x L, eq. (20) structure

% Explicitly diagnose numerical rank without issuing MATLAB's backslash
% rank-deficiency warning. Z_S loses rank when two paths are assigned the
% same (EM state, grid angle) pair by Steps 3-4, which is exactly what
% happens in the T=4 breakdown of Fig. (c), where the EM coherence
% of eq. (16) reaches 0.985.
svals=svd(ZS,'econ');
if isempty(svals)
    tol_rank=0;
    rank_ZS=0;
    cond_ZS=Inf;
else
    tol_rank=max(size(ZS))*eps(max(svals));
    rank_ZS=sum(svals>tol_rank);
    if svals(end)>0
        cond_ZS=svals(1)/svals(end);
    else
        cond_ZS=Inf;
    end
end

% Minimum-norm LS remains well defined if recovered atoms collide, so a
% rank-deficient realization degrades gracefully instead of returning
% Inf/NaN and poisoning the Monte-Carlo average. lsqminnorm needs
% R2017b+; pinv is the fallback and is equivalent at these sizes.
if exist('lsqminnorm','file')==2
    B_hat=(lsqminnorm(ZS,Y3.',tol_rank)).';
else
    B_hat=(pinv(ZS,tol_rank)*Y3.').';
end

%% Step 7: path-gain recovery, eq. (27)
% B = diag(s) C with diag(s) invertible whenever s_p != 0, so the
% structured LS problem (27) has the exact closed-form solution
% C = diag(s)^-1 B: divide each row of B_hat by its pilot symbol.
if any(abs(data.s)<=eps)
    error('The single-user pilot sequence contains a zero entry.');
end
alpha_hat=B_hat./data.s;

%% Outputs
phys=struct();
phys.HEM_hat=HEM_hat;
phys.A_hat=A_hat;
phys.theta_hat=theta_hat;
phys.TH_phys=TH_phys;
phys.TA_phys=TA_phys;
phys.B_hat=B_hat;
phys.alpha_hat=alpha_hat;

% Diagnostics for difficult realizations.
phys.rank_ZS=rank_ZS;
phys.rank_deficient_ZS=(rank_ZS<R);
phys.cond_ZS=cond_ZS;
phys.singular_values_ZS=svals;
end

function phys = trace_somp_baseline(data,base,R)
%TRACE_SOMP_BASELINE Greedy joint EM/AoA baseline (Section V).
%
%   phys = TRACE_SOMP_BASELINE(data,base,R)
%
% Simultaneous orthogonal matching pursuit over the joint EM/angular
% dictionary built by TRACE_SOMP_PREPARE. This is the dictionary-based
% competitor the paper compares against: it represents the classical way
% of posing hybrid-array channel estimation, extended to the tri-hybrid
% model, and it treats the EM and spatial components as one flat sparse
% recovery instead of exploiting the multilinear coupling.
%
% ALGORITHM
%   Observation matrix V = Y_(3)^T in C^((NRF*T) x P): the P pilot slots
%   are the multiple measurement vectors, and they share one support
%   because the AoAs and EM signatures are constant over the slots while
%   only the gains vary. For r = 1,...,R:
%     1. score every atom by its total correlation energy over the P
%        columns of the residual, sum_p |d^H res_p|^2 (the "simultaneous"
%        part of SOMP)
%     2. exclude atoms already chosen, take the maximizer
%     3. re-solve the LS problem on the enlarged support and update the
%        residual
%   A final LS refinement on the selected support gives the coefficients.
%
% WHAT THE COMPARISON SHOWS (Fig. (b))
%   SOMP is the more accurate estimator below the crossover at about
%   7 dB, because a greedy support search is robust when the tensor fit
%   is noise-dominated. Above it, SOMP saturates at a grid-induced error
%   floor -- its accuracy is bounded by how well a grid atom can
%   represent an off-grid AoA -- whereas TRACE keeps the SNR slope of the
%   oracle-support bound. The two methods therefore fail in different
%   ways, which is the point of plotting both.
%
% INPUTS
%   data  realization structure (needs .Y, .s, .par.Nr, .par.M)
%   base  dictionary from TRACE_SOMP_PREPARE
%   R     number of paths to extract (R = L; the true L is assumed known,
%         as it is for TRACE, so neither method is penalized for model
%         order selection)
%
% OUTPUTS
%   phys.A_hat, phys.HEM_hat, phys.B_hat, phys.alpha_hat, phys.support
%   in the same format as TRACE_PHYSICAL_RECOVERY, so both estimators can
%   be scored by the same TRACE_CHANNEL_NMSE call.
%
% PAPER CROSS-REFERENCE
%   Section V, greedy dictionary baseline; Fig. (b)
%
% See also TRACE_SOMP_PREPARE, TRACE_ORACLE_BASELINE, TRACE_CHANNEL_NMSE.

V=trace_unfold(data.Y,3).';                % (NRF*T) x P, eq. (7)/(20)
res=V;
supp=zeros(1,R);

for r=1:R
    % Simultaneous scoring: total correlation energy across all P slots.
    score=sum(abs(base.Dnorm'*res).^2,2);
    if r>1
        score(supp(1:r-1))=-inf;           % no atom selected twice
    end

    [~,supp(r)]=max(score);
    X=base.D(:,supp(1:r))\V;               % LS on the current support
    res=V-base.D(:,supp(1:r))*X;
end

% Final coefficient refinement on selected support.
X=base.D(:,supp)\V;
B_hat=X.';

% Decode joint EM/AoA dictionary indices: column index of D is
% (m-1)*G + g, so g is the fast index and m the slow one.
g=mod(supp-1,base.G)+1;
m=floor((supp-1)/base.G)+1;

A_hat=ula_steering(data.par.Nr,base.theta_grid(g));
HEM_hat=zeros(data.par.M,R);
HEM_hat(sub2ind(size(HEM_hat),m,1:R))=1;

% Single-user pilot removal, as in eq. (27).
alpha_hat=B_hat./data.s;

phys=struct();
phys.A_hat=A_hat;
phys.HEM_hat=HEM_hat;
phys.B_hat=B_hat;
phys.alpha_hat=alpha_hat;
phys.support=supp;
end

function nmse = trace_channel_nmse(A0,H0,alpha0,A1,H1,alpha1)
%TRACE_CHANNEL_NMSE The performance metric reported in Section V.
%
%   nmse = TRACE_CHANNEL_NMSE(A0,H0,alpha0,A1,H1,alpha1)
%
% Computes the channel NMSE defined in Section V,
%
%     NMSE = sum_p || hbar_p - hbar_p_hat ||_2^2 / sum_p || hbar_p ||_2^2,
%
% where hbar_p is the extended (EM x spatial) channel of eq. (3),
%
%     hbar_p = (A(theta) <> H_EMbar) alpha_p
%            = sum_l alpha_{p,l} (a(theta_l) kron h_{EM,l}).
%
% This is the ONLY headline metric in the paper. It is an end-to-end
% measure: it scores the reconstructed physical channel, not the latent
% CP factors, so it is unaffected by how well any individual factor was
% recovered as long as their product is right.
%
% PERMUTATION INVARIANCE
%   The sum over paths is invariant to a COMMON permutation of the
%   columns of (A, H_EMbar, alpha), which is exactly the CP ambiguity.
%   No path matching, and in particular no oracle knowledge, is needed to
%   evaluate it. This is why the metric can be computed directly on the
%   estimator output, unlike the per-factor NMSEs in TRACE_EVAL.
%
% RATIO OF SUMS, NOT MEAN OF RATIOS
%   Within one realization the P slots are accumulated in the numerator
%   and denominator separately. Across realizations, the reproduction
%   scripts store the per-realization value and the paper reports the
%   MEDIAN over the MC = 200 draws; the median is used because the
%   distribution is heavy-tailed near the low-SNR threshold, where a few
%   failed realizations (NMSE > 1) would otherwise dominate a mean.
%
% INPUTS
%   A0,H0,alpha0   true      A(theta) [Nr x L], H_EMbar [M x L], C [P x L]
%   A1,H1,alpha1   estimated counterparts, same sizes
%
% PAPER CROSS-REFERENCE
%   eq. (3) for hbar_p; Section V, "Performance is measured by the
%   channel NMSE ..."
%
% See also TRACE_EVAL, TRACE_PHYSICAL_RECOVERY.

G0=khatrirao(A0,H0);
G1=khatrirao(A1,H1);
P=size(alpha0,1);

num=0;
den=0;

for p=1:P
    htrue=G0*alpha0(p,:).';
    hhat =G1*alpha1(p,:).';

    num=num+norm(htrue-hhat)^2;
    den=den+norm(htrue)^2;
end

nmse=num/max(den,eps);
end

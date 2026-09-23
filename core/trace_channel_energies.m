function [num,den] = trace_channel_energies(A0,H0,alpha0,A1,H1,alpha1)
%TRACE_CHANNEL_ENERGIES Squared channel error and true channel energy.
%
%   [num,den] = TRACE_CHANNEL_ENERGIES(A0,H0,alpha0,A1,H1,alpha1)
%
% Same quantity as TRACE_CHANNEL_NMSE, but returning the numerator and
% the denominator separately instead of their ratio. The whole set of P
% channels is formed at once,
%
%     Hbar = (A(theta) <> H_EMbar) C^T,          eq. (3) stacked over p,
%
% so num = ||Hbar - Hbar_hat||_F^2 and den = ||Hbar||_F^2.
%
% WHY SEPARATE ENERGIES
%   Keeping them apart lets one run support both statistics used in the
%   paper: the per-trial NMSE num/den, whose MEDIAN over realizations is
%   what the figures plot and from which the NMSE>1 failure rates are
%   computed; and the classical accumulated Monte-Carlo NMSE
%   sum(num)/sum(den). Storing only the ratio would make the latter
%   unrecoverable.
%
% This function is used by REPRO_FIG_C_PROBING_STATES.
%
% PAPER CROSS-REFERENCE
%   eq. (3); Section V, channel NMSE definition
%
% See also TRACE_CHANNEL_NMSE.

G0 = khatrirao(A0,H0);
G1 = khatrirao(A1,H1);

Htrue = G0 * alpha0.';
Hhat  = G1 * alpha1.';

num = norm(Htrue-Hhat,'fro')^2;
den = norm(Htrue,'fro')^2;
end

function phys = trace_oracle_baseline(data)
%TRACE_ORACLE_BASELINE Oracle-support LS bound (Section V).
%
%   phys = TRACE_ORACLE_BASELINE(data)
%
% A REFERENCE BOUND, NOT A COMPETING ESTIMATOR. The true A(theta) and
% H_EMbar are handed to the estimator, so the AoAs and EM signatures cost
% nothing; only the pilot-weighted gain matrix B (eq. (10)) is estimated,
% by least squares from the mode-3 unfolding:
%
%     B_hat = Y_(3) [(T_H <> T_A)^T]^dagger,
%
% i.e. eq. (26) with the TRUE structured factors in place of the
% recovered ones, followed by C = diag(s)^-1 B as in eq. (27).
%
% WHAT IT ISOLATES
%   Everything the oracle curve leaves is the irreducible cost of
%   estimating L*P complex gains from NRF*T*P noisy observations. The gap
%   between TRACE and this curve is therefore exactly the price of not
%   knowing the support -- it is the quantity Section V quotes as the
%   "TRACE/oracle offset", and the fact that it settles to a small
%   constant factor at high SNR is the paper's central accuracy claim.
%
%   The oracle also inherits the 2x-per-3-dB slope of an unbiased LS
%   estimator, which is why a curve that tracks it in slope (TRACE) and
%   one that does not (SOMP, which hits a grid floor) are so easy to
%   distinguish in Fig. (b).
%
% INPUT
%   data   realization structure; uses .Y, .TA, .TH (TRUE factors),
%          .A, .HEM (true physical parameters) and .s
%
% OUTPUT
%   phys   same fields as TRACE_PHYSICAL_RECOVERY, so the same
%          TRACE_CHANNEL_NMSE call scores it
%
% PAPER CROSS-REFERENCE
%   eqs. (26), (27); Section V, oracle-support bound
%
% See also TRACE_PHYSICAL_RECOVERY, TRACE_SOMP_BASELINE.

Y3=trace_unfold(data.Y,3);
ZS=khatrirao(data.TH,data.TA);             % TRUE Z_S, eq. (20)

% lsqminnorm (R2017b+) with pinv as the fallback; identical solution here
% because Z_S is full rank whenever the true angles and EM states are
% distinct, which they are by construction.
if exist('lsqminnorm','file')==2
    B_hat=(lsqminnorm(ZS,Y3.')).';
else
    B_hat=(pinv(ZS)*Y3.').';
end

alpha_hat=B_hat./data.s;                   % eq. (27)

phys=struct();
phys.A_hat=data.A;
phys.HEM_hat=data.HEM;
phys.B_hat=B_hat;
phys.alpha_hat=alpha_hat;
end

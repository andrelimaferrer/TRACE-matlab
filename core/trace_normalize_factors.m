function [TA,TH,TS] = trace_normalize_factors(TA,TH,TS)
%TRACE_NORMALIZE_FACTORS Fix the CP scaling ambiguity without changing the tensor.
%
%   [TA,TH,TS] = TRACE_NORMALIZE_FACTORS(TA,TH,TS)
%
% A PARAFAC model is invariant under reciprocal column scalings,
%
%   t_{A,l} o t_{H,l} o t_{S,l}
%     = (t_{A,l}/c) o (t_{H,l}/d) o (c*d*t_{S,l}),
%
% which is one of the two ambiguities discussed in Section V of the paper
% (the other being column permutation). Left unchecked during ALS, the
% norms of T_A and T_H can drift over many orders of magnitude while the
% fit stays constant, which destroys the conditioning of the backslash
% solves in eq. (21). This function removes the drift by giving T_A and
% T_H unit-norm columns and absorbing all the scale into T_S.
%
% The reconstructed tensor is bit-for-bit unchanged up to floating-point
% rounding: this is a re-parameterization, not an estimation step.
%
% The scaling ambiguity is immaterial for the reported results because
% the physical recovery stage (eq. (25)-(26)) re-references the third
% factor to the known training matrices, and the channel NMSE of
% Section V is computed on the reconstructed channel, not on the factors.
%
% PAPER CROSS-REFERENCE
%   Scaling ambiguity: Section V, first paragraph; eqs. (11), (15), (26)
%
% See also TRACE_ALS, TRACE_PHYSICAL_RECOVERY.

nA=max(vecnorm(TA),eps); TA=TA./nA; TS=TS.*nA;
nH=max(vecnorm(TH),eps); TH=TH./nH; TS=TS.*nH;
end

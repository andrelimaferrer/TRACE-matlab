function init = trace_random_init(Y,R)
%TRACE_RANDOM_INIT Random starting point for the ALS cycle.
%
%   init = TRACE_RANDOM_INIT(Y,R)
%
% Produces one independent random initialization of the three factors.
% The problem (17) is non-convex, so Step 2 of Algorithm 1 restarts the
% ALS cycle from N_s - 1 such points (N_s = 5 in the paper) in addition
% to the SVD start, and keeps the smallest-residual fit.
%
% The factors are drawn i.i.d. complex Gaussian, normalized by
% TRACE_NORMALIZE_FACTORS, and then the whole model is rescaled so that
%
%     || [[T_A,T_H,T_S]] ||_F  =  || Y ||_F .
%
% Matching the initial energy to the data matters: an initialization
% whose norm is orders of magnitude off the observation spends the first
% ALS sweeps simply correcting the overall scale, and with a finite
% iteration budget (I_max = 200) that can cost accuracy.
%
% INPUTS
%   Y    NRF x T x P pilot tensor
%   R    CP rank; R = L for the single-user model
%
% OUTPUT
%   init.TA, init.TH, init.TS, init.method
%
% PAPER CROSS-REFERENCE
%   Algorithm 1, Step 2 (multi-start); eq. (17)
%
% See also TRACE_ALS, TRACE_NORMALIZE_FACTORS.

[NRF,T,P]=size(Y);
TA=randn(NRF,R)+1j*randn(NRF,R);
TH=randn(T,R)+1j*randn(T,R);
TS=randn(P,R)+1j*randn(P,R);
[TA,TH,TS]=trace_normalize_factors(TA,TH,TS);
% Match the overall initial tensor norm to the observation norm.
Y0=trace_cpdgen(TA,TH,TS);
TS=TS*(norm(Y(:))/max(norm(Y0(:)),eps));
init=struct('TA',TA,'TH',TH,'TS',TS,'method','random');
end

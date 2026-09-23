function Y = trace_cpdgen(TA,TH,TS)
%TRACE_CPDGEN Build the TRACE pilot tensor from its three PARAFAC factors.
%
%   Y = TRACE_CPDGEN(TA,TH,TS)
%
% Evaluates the noiseless part of eq. (11),
%
%   Y = I_{3,L} x_1 T_A x_2 T_H x_3 T_S,
%
% equivalently the sum of L rank-one outer products of eq. (15),
%
%   Y = sum_{l=1}^{L} t_{A,l} o t_{H,l} o t_{S,l},
%
% with one rank-one term per propagation path.
%
% INPUTS
%   TA   NRF x R,  RF-projected steering factor  (eq. (12))
%   TH   T   x R,  EM-domain angular factor      (eq. (13))
%   TS   P   x R,  pilot-weighted fading factor  (eq. (14))
%
% OUTPUT
%   Y    NRF x T x P
%
% IMPLEMENTATION
%   Built from the mode-1 matricization Y_(1) = T_A (T_S <> T_H)^T of
%   eq. (18) and reshaped, which is far faster than an explicit loop over
%   the R rank-one terms. The reshape ordering is the inverse of the one
%   in TRACE_UNFOLD(.,1), so the two functions are exact inverses.
%
% PAPER CROSS-REFERENCE
%   eqs. (11), (15), (18)
%
% See also TRACE_UNFOLD, KHATRIRAO, TRACE_GENERATE_CLEAN_DATA.

R=size(TA,2);
if size(TH,2)~=R || size(TS,2)~=R, error('All factors must have R columns.'); end
Y1=TA*khatrirao(TS,TH).';                 % NRF x (T*P)
Y=reshape(Y1,size(TA,1),size(TH,1),size(TS,1));
end

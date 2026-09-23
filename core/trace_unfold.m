function Ymat = trace_unfold(Y,mode)
%TRACE_UNFOLD Mode-n unfolding of the TRACE pilot tensor, paper convention.
%
%   Ymat = TRACE_UNFOLD(Y,mode)   with   Y in C^(NRF x T x P)
%
% Returns the matricizations Y_(1), Y_(2), Y_(3) used in eqs. (18)-(20):
%
%   mode 1 : Y_(1) = T_A Z_A^T,   Z_A = T_S <> T_H      NRF x (T*P)
%   mode 2 : Y_(2) = T_H Z_H^T,   Z_H = T_S <> T_A      T   x (NRF*P)
%   mode 3 : Y_(3) = T_S Z_S^T,   Z_S = T_H <> T_A      P   x (NRF*T)
%
% ORDERING (this is the part that must be got right)
%   MATLAB is column-major, so reshape(Y,sz(1),[]) lists the second index
%   fastest and the third slowest. Combined with the khatrirao convention
%   (first argument = slow), this gives exactly the pairings above. Mode 2
%   and mode 3 first permute the target mode to the front, so the
%   remaining two modes keep their natural (fast,slow) order:
%
%       mode 2 : permute [2 1 3]  ->  columns ordered (NRF fast, P slow)
%       mode 3 : permute [3 1 2]  ->  columns ordered (NRF fast, T slow)
%
%   Note that mode 3 reproduces the stacking of eq. (7): the RF-chain
%   index is fast and the EM probing index t is slow, so Y_(3)^T is the
%   matrix Y of eq. (7) and eq. (8) is its factorization.
%
% INPUTS
%   Y     NRF x T x P pilot tensor (eq. (11))
%   mode  1, 2 or 3
%
% PAPER CROSS-REFERENCE
%   eqs. (7), (8), (18), (19), (20), (26)
%
% See also KHATRIRAO, TRACE_CPDGEN, TRACE_ALS.

sz=size(Y);
switch mode
    case 1, Ymat=reshape(Y,sz(1),[]);                    % NRF x TP
    case 2, Ymat=reshape(permute(Y,[2 1 3]),sz(2),[]);  % T x NRF*P
    case 3, Ymat=reshape(permute(Y,[3 1 2]),sz(3),[]);  % P x NRF*T
    otherwise, error('mode must be 1, 2, or 3.');
end
end

function KR = khatrirao(A,B)
%KHATRIRAO Column-wise Kronecker (Khatri-Rao) product, written "<>" in the paper.
%
%   KR = KHATRIRAO(A,B)   with   KR(:,r) = kron(A(:,r),B(:,r))
%
% A is I x R, B is J x R, and KR is (I*J) x R.
%
% CONVENTION WARNING
%   The paper writes the Khatri-Rao product of the EM and steering factors
%   as (F_EMbar H_EMbar) <> (F_RF A(theta)) in eq. (8), where the LEFT
%   operand indexes the SLOW dimension. This function follows the MATLAB
%   convention kron(A(:,r),B(:,r)), in which the FIRST argument is also
%   the slow one. Hence the paper expression
%
%       Z_S = T_H <> T_A            (eq. (20))
%
%   is coded as  khatrirao(TH,TA), and similarly
%
%       Z_A = T_S <> T_H            (eq. (18))  ->  khatrirao(TS,TH)
%       Z_H = T_S <> T_A            (eq. (19))  ->  khatrirao(TS,TA)
%
%   The ordering matters: it must match the unfolding convention used by
%   TRACE_UNFOLD, otherwise the ALS updates fit a permuted tensor.
%
% IMPLEMENTATION
%   Fully vectorized via implicit expansion; no loop over the R columns.
%   The reshape produces the elements of B varying fastest, which is the
%   MATLAB column-major layout of kron(A(:,r),B(:,r)).
%
% PAPER CROSS-REFERENCE
%   eqs. (3), (8), (11), (18)-(20), (26)
%
% See also TRACE_UNFOLD, TRACE_CPDGEN, TRACE_ALS.

if size(A,2)~=size(B,2), error('A and B must have the same number of columns.'); end
[m,R]=size(A); n=size(B,1);
KR=reshape(reshape(B,n,1,R).*reshape(A,1,m,R),n*m,R);
end

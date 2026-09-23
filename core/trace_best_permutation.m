function perm = trace_best_permutation(Score)
%TRACE_BEST_PERMUTATION Resolve the CP column permutation for diagnostics.
%
%   perm = TRACE_BEST_PERMUTATION(Score)
%
% Finds the permutation that maximizes sum_r Score(r,perm(r)), i.e. the
% best one-to-one assignment between true and estimated CP columns.
%
% THIS IS A SIMULATION-ONLY UTILITY
%   The estimator itself never needs it. The headline channel NMSE of
%   Section V is permutation-invariant (see TRACE_CHANNEL_NMSE) and
%   Steps 3-7 of Algorithm 1 are permutation-equivariant. Column matching
%   is required only to report PER-PATH quantities -- the individual
%   factor NMSEs and the EM support error rate in TRACE_EVAL -- which are
%   diagnostics, not results.
%
% THREE IMPLEMENTATIONS, PICKED AT RUN TIME
%   1. matchpairs (R2019a+): exact Hungarian algorithm, O(R^3). The cost
%      is negated because matchpairs minimizes.
%   2. R <= 8: exhaustive search over all R! permutations. Exact, and at
%      R <= 8 (8! = 40320) cheaper than setting up anything smarter. All
%      configurations in the paper have L <= 6, so this branch is what
%      runs on an older MATLAB or on Octave.
%   3. R > 8 with no matchpairs: greedy row-by-row assignment. Not
%      guaranteed optimal; reached only outside the paper's parameter
%      range, and it degrades a diagnostic, never a reported result.
%
% INPUT
%   Score  R x R nonnegative matching score; Score(r,c) is the similarity
%          between true column r and estimated column c
%
% OUTPUT
%   perm   1 x R, estimated column perm(r) matches true column r
%
% PAPER CROSS-REFERENCE
%   Permutation ambiguity: Section VI, first paragraph
%
% See also TRACE_EVAL, TRACE_CHANNEL_NMSE.

R=size(Score,1);
if exist('matchpairs','file')==2
    pairs=matchpairs(-Score,1e6); pairs=sortrows(pairs,1); perm=pairs(:,2).'; return;
end
if R<=8
    Pm=perms(1:R); vals=zeros(size(Pm,1),1);
    rows=1:R;
    for q=1:size(Pm,1), vals(q)=sum(Score(sub2ind([R R],rows,Pm(q,:)))); end
    [~,q]=max(vals); perm=Pm(q,:); return;
end
% Greedy fallback for large R when matchpairs is unavailable.
perm=zeros(1,R); used=false(1,R);
for r=1:R
    s=Score(r,:); s(used)=-inf; [~,perm(r)]=max(s); used(perm(r))=true;
end
end

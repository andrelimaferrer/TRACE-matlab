function setup_trace()
%SETUP_TRACE Put the TRACE package on the MATLAB/Octave path.
%
%   setup_trace
%
% Run once per session, from anywhere. The function locates the package
% from its own file location, so it does not matter what the current
% working directory is.
%
% DEPENDENCIES
%   None. Everything ships with this package and runs on a base MATLAB
%   installation (no toolboxes) and on GNU Octave. In particular:
%
%     * The CP/PARAFAC solver used in the paper (trace_als) is implemented
%       here from scratch -- no Tensorlab and no Tensor Toolbox.
%     * khatrirao, trace_unfold and trace_cpdgen replace the corresponding
%       third-party helpers.
%     * Optional functions are detected at run time and replaced by
%       equivalent fallbacks when absent:
%         lsqminnorm     (R2017b+) -> pinv
%         matchpairs     (R2019a+) -> exhaustive/greedy matching
%         exportgraphics (R2020a+) -> print
%       The fallbacks are mathematically equivalent for the problem sizes
%       used in the paper; only speed differs. Tensorlab's cpderr is
%       replaced by the factor alignment built into trace_eval, and
%       dftmtx (Signal Processing Toolbox) by a local helper.
%
% REFERENCE
%   A. L. F. de Almeida et al., "TRACE: A Trilinear Channel Estimation
%   Framework for Tri-Hybrid Architecture."
%   Equation numbers quoted in the source headers refer to that paper.
%
% See also DEMO_TRACE, TRACE_ALS, TRACE_PHYSICAL_RECOVERY.

root = fileparts(mfilename('fullpath'));

addpath(root);
addpath(fullfile(root,'core'));
addpath(fullfile(root,'baselines'));
addpath(fullfile(root,'reproduce'));

% Sanity check: the package must be self-consistent.
required = {'trace_als','trace_generate_clean_data','trace_physical_recovery', ...
            'trace_design_training','khatrirao','trace_unfold'};
missing = required(cellfun(@(f) exist(f,'file')~=2, required));
if ~isempty(missing)
    error('TRACE:IncompleteInstall', ...
        ['The following package functions were not found: %s\n' ...
         'Make sure core/ and baselines/ were not removed from %s.'], ...
        strjoin(missing,', '), root);
end

fprintf('TRACE single-user package ready (%s).\n', root);
fprintf('Solver: TRACE-ALS (SVD init + QR/backslash, multi-start).\n');
fprintf('Start with:  demo_trace          -- one realization, ~2 s\n');
fprintf('             help reproduce      -- paper figure scripts\n');

end

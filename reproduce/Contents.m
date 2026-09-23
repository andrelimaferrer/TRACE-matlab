% REPRODUCE  Scripts that regenerate the paper's results figure.
%
% Run any of the three simulation scripts, then the renderer. Each script
% fixes its own seed and writes one file into results/.
%
%   repro_fig_a_training_dimensions  - Panel (a): NMSE vs SNR for three
%                                      training dimension choices,
%                                      non-compressive. ~10-25 min.
%   repro_fig_b_baselines            - Panel (b): TRACE vs greedy SOMP and
%                                      the oracle bound, compressive.
%                                      ~30-60 min.
%   repro_fig_c_probing_states       - Panel (c): NMSE vs the number of EM
%                                      probing states T at 10 dB. ~20-40 min.
%   make_paper_figure                - Draws the three-panel figure from
%                                      results/. Seconds.
%
% Times are for the published MC = 200. Set MC = 20 at the top of a script
% for a quick look.
%
% The shipped results/ already contains the published data, so
% make_paper_figure works without running any simulation.
%
% See also SETUP_TRACE, DEMO_TRACE.

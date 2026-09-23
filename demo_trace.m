%% DEMO_TRACE  Minimal end-to-end run of the TRACE estimator.
%
% One channel realization, one SNR point, the full Algorithm 1 pipeline.
% Runs in a couple of seconds and is the right place to start reading the
% code. Every step below is annotated with the equation it implements.
%
% Reference:
%   A. L. F. de Almeida et al., "TRACE: A Trilinear Channel Estimation
%   Framework for Tri-Hybrid Architecture."
%
% To reproduce the published figures instead, see reproduce/.

clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir); setup_trace;

rng(7);                            % reproducible run

%% ---------------------------------------------------------------------
% System configuration. K = 1 is implicit throughout and R = L (eq. (15)).
%% ---------------------------------------------------------------------
par.L   = 4;                       % propagation paths = tensor rank
par.Nr  = 4;                       % receive antennas
par.NRF = 4;                       % RF chains
par.M   = 16;                      % EM reconfiguration states
par.T   = 24;                      % EM probing states per slot
par.P   = 8;                       % pilot slots
par.FEMbar_type = 'random-phase';  % eq. (9); see TRACE_DESIGN_TRAINING
par.FRF_type    = 'random-phase';
par.min_angle_sep_deg = 15;

SNRdB = 15;
R = par.L;

% Check the design rule of eq. (29) before spending any time on it.
kH = min(par.T,par.L); kA = min(par.NRF,par.L); kS = min(par.P,par.L);
relsym = '<'; if kH+kA+kS >= 2*par.L+2, relsym = '>='; end
fprintf('Kruskal check (eq. 29): %d + %d + %d = %d %s 2L+2 = %d\n', ...
    kH,kA,kS,kH+kA+kS,relsym,2*par.L+2);

%% ---------------------------------------------------------------------
% Training design and one noisy realization
%% ---------------------------------------------------------------------
[training,trainDiag] = trace_design_training(par);       % eqs. (9), (16)
fprintf('Training: mu_EM = %.3f (eq. 16), duplicate EM pairs = %d, min RF gain = %.3f\n', ...
    trainDiag.mu_EM,trainDiag.duplicate_EM_pairs,trainDiag.min_RF_gain_rel);

clean = trace_generate_clean_data(par,training);         % eqs. (2)-(3), (11)-(14)
data  = trace_add_noise(clean,SNRdB);                    % eq. (11), noise term

%% ---------------------------------------------------------------------
% Algorithm 1, Steps 1-2: PARAFAC fit by multi-start ALS
%% ---------------------------------------------------------------------
oa = struct();
oa.init_method = 'svd';            % Step 1
oa.nstarts     = 5;                % N_s = 5, as in the paper
oa.maxit       = 200;              % I_max
oa.tol         = 1e-8;             % epsilon
oa.verbose     = true;

est = trace_als(data.Y,R,oa);                            % eqs. (17)-(21)

%% ---------------------------------------------------------------------
% Algorithm 1, Steps 3-7: physical recovery
%% ---------------------------------------------------------------------
thetaGrid = linspace(-pi/2,pi/2,5001);                   % the grid Theta
rec  = trace_prepare_recovery(training,par.Nr,thetaGrid);
phys = trace_physical_recovery(est,data,rec);            % eqs. (22)-(27)

%% ---------------------------------------------------------------------
% Scoring (needs ground truth, so simulation only)
%% ---------------------------------------------------------------------
met = trace_eval(data,est,phys);

fprintf('\nTRACE single-user validation at SNR = %g dB\n',SNRdB);
fprintf('  Tensor fitting NMSE (eq. 17) : %.4e\n',est.tensor_nmse);
fprintf('  T_A factor NMSE     (eq. 12) : %.4e\n',met.nmse_TA);
fprintf('  T_H factor NMSE     (eq. 13) : %.4e\n',met.nmse_TH);
fprintf('  T_S factor NMSE     (eq. 14) : %.4e\n',met.nmse_TS);
fprintf('  H_EM NMSE           (eq.  3) : %.4e\n',met.nmse_HEM);
fprintf('  A(theta) NMSE       (eq.  2) : %.4e\n',met.nmse_A);
fprintf('  path gain NMSE      (eq. 27) : %.4e\n',met.nmse_alpha);
fprintf('  EM support error rate        : %.2f %%\n',100*met.HEM_support_error_rate);
fprintf('  rank(Z_S)           (eq. 26) : %d / %d\n',met.rank_ZS,R);
fprintf('  CHANNEL NMSE (Section V)     : %.4e\n',met.channel_nmse);

% AoA accuracy, matched to the true paths for readability only.
theta_true = sort(data.theta)*180/pi;
theta_est  = sort(phys.theta_hat)*180/pi;
fprintf('\n  True AoAs (deg) : %s\n',num2str(theta_true,'%8.2f'));
fprintf('  Est. AoAs (deg) : %s\n',num2str(theta_est,'%8.2f'));

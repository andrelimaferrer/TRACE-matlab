%% REPRO_FIG_A_TRAINING_DIMENSIONS
% Fig. 2(a) of the paper: channel NMSE versus SNR for three training
% dimension choices, in a NON-COMPRESSIVE setting (N_RF = N_r, T >= M).
%
% WHAT THE PANEL SHOWS
%   Pilot slots beyond the tensor rank remove the low-SNR threshold effect
%   seen at P = L, while extra EM probing slots give a uniform gain. The
%   P = L = 4 curve is the interesting one: eq. (29) is satisfied with
%   zero margin there (min(T,L)+min(NRF,L)+min(P,L) = 4+4+4 = 12 >= 10),
%   yet the estimate breaks down at low SNR -- which is Remark 3 of the
%   paper: identifiability is not conditioning.
%
% CONFIGURATION (Section V)
%   N_r = N_RF = 4, M = 16, L = 4, minimum AoA separation 15 degrees
%   Cfg rows [P, T, N_RF]:  (4,16,4), (8,16,4), (8,24,4)
%   SNR 0:3:21 dB, MC = 200 realizations, N_s = 5 ALS starts
%
% RUNTIME
%   Roughly 10-25 min on a current laptop for MC = 200. Set MC = 20 for a
%   quick look; the curve shapes are already clear, only the tails move.
%
% OUTPUT
%   results/figure_a_data.mat with the per-realization NMSE matrices, so
%   the medians and the NMSE>1 failure rates quoted in Section V can be
%   recomputed without re-running.
%
% See also MAKE_PAPER_FIGURE, TRACE_ALS, TRACE_PHYSICAL_RECOVERY.

clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
addpath(fileparts(thisDir)); setup_trace;

rng(11);                           % fixed seed for the whole experiment

%% ---------------------------------------------------------------------
% Fixed parameters
%% ---------------------------------------------------------------------
base.L   = 4;
base.Nr  = 4;
base.M   = 16;
base.FEMbar_type = 'random-phase';
base.FRF_type    = 'random-phase';
base.min_angle_sep_deg = 15;

SNRdB   = 0:3:21;
MC      = 200;                     % Monte-Carlo realizations
NSTARTS = 5;                       % N_s in Algorithm 1, Step 2

% Rows: [P, T, NRF]
Cfg = [4 16 4;
       8 16 4;
       8 24 4];

labels = {'P=L=4, T=16, N_{RF}=4', ...
          'P=8, T=16, N_{RF}=4', ...
          'P=8, T=24, N_{RF}=4'};

all_nmse = cell(size(Cfg,1),1);    % per-realization NMSE, one cell per config
curves   = zeros(numel(SNRdB),size(Cfg,1));   % medians, for a quick plot

%% ---------------------------------------------------------------------
% Monte-Carlo sweep
%% ---------------------------------------------------------------------
for c = 1:size(Cfg,1)

    par = base;
    par.P   = Cfg(c,1);
    par.T   = Cfg(c,2);
    par.NRF = Cfg(c,3);
    R = par.L;

    % Training is drawn ONCE per configuration and reused for every
    % realization and SNR point: a fixed training protocol (Section V).
    [training,trainDiag] = trace_design_training(par);
    fprintf('Cfg %d: mu_EM=%.3f, duplicate EM pairs=%d, min RF gain=%.3f\n', ...
        c,trainDiag.mu_EM,trainDiag.duplicate_EM_pairs,trainDiag.min_RF_gain_rel);

    thetaGrid = linspace(-pi/2,pi/2,2001);
    rec = trace_prepare_recovery(training,par.Nr,thetaGrid);

    tmp = zeros(numel(SNRdB),MC);

    for mc = 1:MC

        % One channel realization, reused across all SNR points so the
        % curves are paired.
        clean = trace_generate_clean_data(par,training);

        for is = 1:numel(SNRdB)
            data = trace_add_noise(clean,SNRdB(is));

            oa = struct();
            oa.init_method = 'svd';          % Algorithm 1, Step 1
            oa.nstarts     = NSTARTS;
            oa.maxit       = 200;
            oa.tol         = 1e-8;
            oa.verbose     = false;

            est  = trace_als(data.Y,R,oa);              % eqs. (17)-(21)
            phys = trace_physical_recovery(est,data,rec); % eqs. (22)-(27)
            met  = trace_eval(data,est,phys);

            tmp(is,mc) = met.channel_nmse;
        end

        if mod(mc,25)==0
            fprintf('  Configuration %d/%d, MC %d/%d\n',c,size(Cfg,1),mc,MC);
        end
    end

    all_nmse{c} = tmp;
    curves(:,c) = median(tmp,2);   % the paper reports the MEDIAN
end

%% ---------------------------------------------------------------------
% Quick look (the publication figure is drawn by make_paper_figure)
%% ---------------------------------------------------------------------
figure;
semilogy(SNRdB,curves(:,1),'-o','LineWidth',1.6,'MarkerSize',7); hold on;
semilogy(SNRdB,curves(:,2),'-s','LineWidth',1.6,'MarkerSize',7);
semilogy(SNRdB,curves(:,3),'-d','LineWidth',1.6,'MarkerSize',7);
grid on; box on;
xlabel('SNR (dB)'); ylabel('NMSE');
legend(labels,'Location','southwest');
title('TRACE: impact of training dimensions');
set(gca,'FontSize',11);

%% ---------------------------------------------------------------------
% Save
%% ---------------------------------------------------------------------
outDir = fullfile(fileparts(thisDir),'results');
if ~exist(outDir,'dir'), mkdir(outDir); end
save(fullfile(outDir,'figure_a_data.mat'), ...
     'SNRdB','Cfg','labels','base','curves','all_nmse','MC','NSTARTS');

fprintf('\nSaved %s\n',fullfile(outDir,'figure_a_data.mat'));
fprintf('Median NMSE, P=L=4 : %.3e at %g dB -> %.3e at %g dB\n', ...
    curves(1,1),SNRdB(1),curves(end,1),SNRdB(end));
fprintf('NMSE > 1 rate, P=L=4 at %g dB : %.1f %%\n', ...
    SNRdB(1),100*mean(all_nmse{1}(1,:)>1));

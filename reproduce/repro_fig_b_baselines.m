%% REPRO_FIG_B_BASELINES
% Fig. 2(b) of the paper: TRACE against the greedy dictionary baseline
% and the oracle-support bound, in a COMPRESSIVE configuration.
%
% WHAT THE PANEL SHOWS
%   The configuration is genuinely compressive in both reduced domains:
%   N_RF = 3 < N_r = 8 (fewer RF chains than antennas) and T = 10 < M = 24
%   (fewer probing states than reconfiguration states). The greedy
%   baseline is more accurate below the crossover at about 7 dB, but
%   saturates at a grid-induced error floor, whereas TRACE keeps the SNR
%   slope of the oracle-support bound (a factor 2 per 3 dB). The two
%   estimators therefore fail in different ways, which is the point.
%
% CONFIGURATION (Section V)
%   L = 4, N_r = 8, N_RF = 3, M = 24, T = 10, P = 6
%   minimum AoA separation 5 degrees, SNR -5:3:22 dB
%   MC = 200 realizations, N_s = 5 ALS starts
%   Angular grid: 1001 points, THE SAME grid for TRACE and SOMP so that
%   neither method benefits from a finer dictionary.
%
%   Identifiability check, eq. (29):
%     min(T,L)+min(N_RF,L)+min(P,L) = 4+3+4 = 11 >= 2L+2 = 10.  OK
%
% NOTE ON THE PLOTTED RANGE
%   The paper plots this panel from 1 dB upward. The sweep starts at
%   -5 dB so the threshold region is recorded in the saved data, but
%   below about 1 dB all three estimators are in breakdown and the curves
%   carry no information.
%
% RUNTIME
%   Roughly 30-60 min at MC = 200; the SOMP dictionary has M*G = 24024
%   atoms and dominates. Set MC = 20 for a quick look.
%
% OUTPUT
%   results/figure_b_data.mat with the per-realization NMSE matrices.
%
% See also MAKE_PAPER_FIGURE, TRACE_SOMP_BASELINE, TRACE_ORACLE_BASELINE.

clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
addpath(fileparts(thisDir)); setup_trace;

rng(22);

%% ---------------------------------------------------------------------
% Configuration
%% ---------------------------------------------------------------------
par.L   = 4;
par.Nr  = 8;
par.NRF = 3;                       % compressive: NRF < Nr
par.M   = 24;
par.T   = 10;                      % compressive: T < M
par.P   = 6;
par.FEMbar_type = 'random-phase';
par.FRF_type    = 'random-phase';
par.min_angle_sep_deg = 5;

R       = par.L;
SNRdB   = -5:3:22;
MC      = 200;
NSTARTS = 5;

%% ---------------------------------------------------------------------
% Fixed training, shared dictionaries
%% ---------------------------------------------------------------------
[training,trainDiag] = trace_design_training(par);
fprintf('Training: mu_EM=%.3f (eq. 16), duplicate EM pairs=%d, min RF gain=%.3f\n', ...
    trainDiag.mu_EM,trainDiag.duplicate_EM_pairs,trainDiag.min_RF_gain_rel);

thetaGrid  = linspace(-pi/2,pi/2,1001);   % same grid for TRACE and SOMP
rec        = trace_prepare_recovery(training,par.Nr,thetaGrid);
sompCache  = trace_somp_prepare(training,par.Nr,thetaGrid);

N = numel(SNRdB);
als    = zeros(N,MC);
somp   = zeros(N,MC);
oracle = zeros(N,MC);
als_tensor_nmse = zeros(N,MC);
rankdef_als     = false(N,MC);

%% ---------------------------------------------------------------------
% Monte-Carlo sweep
%% ---------------------------------------------------------------------
for mc = 1:MC

    clean = trace_generate_clean_data(par,training);

    for is = 1:N
        data = trace_add_noise(clean,SNRdB(is));

        %% TRACE: ALS (Algorithm 1 Steps 1-2) + recovery (Steps 3-7)
        oa = struct();
        oa.init_method = 'svd';
        oa.nstarts     = NSTARTS;
        oa.maxit       = 200;
        oa.tol         = 1e-8;
        oa.verbose     = false;

        estALS  = trace_als(data.Y,R,oa);
        physALS = trace_physical_recovery(estALS,data,rec);
        metALS  = trace_eval(data,estALS,physALS);

        als(is,mc)             = metALS.channel_nmse;
        als_tensor_nmse(is,mc) = estALS.tensor_nmse;
        rankdef_als(is,mc)     = metALS.rank_deficient_ZS;

        %% Greedy joint EM/AoA dictionary baseline
        ps = trace_somp_baseline(data,sompCache,R);
        somp(is,mc) = trace_channel_nmse( ...
            data.A,data.HEM,data.alpha, ...
            ps.A_hat,ps.HEM_hat,ps.alpha_hat);

        %% Oracle-support LS reference bound
        po = trace_oracle_baseline(data);
        oracle(is,mc) = trace_channel_nmse( ...
            data.A,data.HEM,data.alpha, ...
            po.A_hat,po.HEM_hat,po.alpha_hat);
    end

    if mod(mc,25)==0
        fprintf('  compressive configuration: MC %d/%d\n',mc,MC);
    end
end

%% ---------------------------------------------------------------------
% Medians (the reported statistic) and a quick look
%% ---------------------------------------------------------------------
med_als    = median(als,2);
med_somp   = median(somp,2);
med_oracle = median(oracle,2);

keep = SNRdB >= 1;                 % the range shown in the paper

figure;
semilogy(SNRdB(keep),med_als(keep),'-o','LineWidth',1.6,'MarkerSize',7); hold on;
semilogy(SNRdB(keep),med_somp(keep),'-.d','LineWidth',1.6,'MarkerSize',7);
semilogy(SNRdB(keep),med_oracle(keep),':^','LineWidth',1.8,'MarkerSize',7);
grid on; box on;
xlabel('SNR (dB)'); ylabel('NMSE');
legend('TRACE-ALS','Joint SOMP','Oracle-support LS','Location','southwest');
title('TRACE versus baselines, compressive configuration');
set(gca,'FontSize',11);

%% ---------------------------------------------------------------------
% Save
%% ---------------------------------------------------------------------
outDir = fullfile(fileparts(thisDir),'results');
if ~exist(outDir,'dir'), mkdir(outDir); end
save(fullfile(outDir,'figure_b_data.mat'), ...
     'par','trainDiag','SNRdB','MC','NSTARTS','R', ...
     'als','somp','oracle','als_tensor_nmse','rankdef_als', ...
     'med_als','med_somp','med_oracle');

fprintf('\nSaved %s\n',fullfile(outDir,'figure_b_data.mat'));

% Reproduce the numbers quoted in Section V for this panel.
j = @(s) find(SNRdB==s,1);
fprintf('SOMP advantage at 1 dB      : %.2fx\n',med_als(j(1))/med_somp(j(1)));
fprintf('crossover ratio at 7 dB     : %.2f (1.00 = crossover)\n', ...
    med_somp(j(7))/med_als(j(7)));
fprintf('TRACE advantage at 22 dB    : %.2fx\n',med_somp(j(22))/med_als(j(22)));
fprintf('TRACE NMSE>1 rate at 1 dB   : %.1f %%\n',100*mean(als(j(1),:)>1));
hi = med_als(SNRdB>=10)./med_oracle(SNRdB>=10);
fprintf('TRACE/oracle offset >=10 dB : %.2f to %.2f\n',min(hi),max(hi));

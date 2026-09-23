%% REPRO_FIG_C_PROBING_STATES
% Fig. 2(c) of the paper: channel NMSE versus the number of EM probing
% states T, at a fixed SNR of 10 dB.
%
% WHAT THE PANEL SHOWS
%   With L = 6, N_RF = 4 and P = 8, the sufficient condition (29)
%
%       min(T,L) + min(N_RF,L) + min(P,L) >= 2L + 2
%       min(T,6) +      4      +     6    >=    14        =>  T >= 4
%
%   is met with ZERO MARGIN at T = 4. It is met, and the estimator still
%   breaks down there: median NMSE 0.82, more than a third of the
%   realizations above 1, an order of magnitude off the oracle. The
%   mechanism is conditioning, not identifiability -- the EM column
%   coherence of eq. (16) is 0.985 at T = 4 against 0.52 at T = 24, so
%   the support detection of eq. (23) can no longer separate the states.
%   Satisfying Kruskal's condition is necessary housekeeping, not a
%   design target.
%
% CONFIGURATION (Section V)
%   L = 6, N_r = 8, N_RF = 4, M = 24, P = 8, minimum AoA separation 5 deg
%   SNR = 10 dB, T in {4, 8, 12, 16, 20, 24}
%   MC = 200 realizations, N_s = 5 ALS starts
%
% THE NESTED-CODEBOOK CONTROL (this is what makes the panel honest)
%   Two things could explain a T-dependence: the number of measurements,
%   or the luck of the particular codebook drawn at each T. To isolate
%   the first, ONE master random-phase EM codebook is drawn at
%   T = max(Tvec) and every smaller design uses its FIRST T ROWS,
%   rescaled by sqrt(T_max/T) so each EM column keeps unit expected norm.
%   The RF combiner is held fixed across all T. Within each Monte-Carlo
%   trial the physical channel (AoAs, EM states, gains, pilots) is also
%   identical across T, so the curves are paired and the only thing that
%   varies is how many probing states were actually used.
%
% RUNTIME
%   Roughly 20-40 min at MC = 200 (6 values of T, four estimators each).
%
% OUTPUT
%   results/figure_c_data.mat, including muEM(T), the coherence of
%   eq. (16) for each codebook, which is the explanatory variable.
%
% NMSE DEFINITION IN THIS SCRIPT
%   The per-trial error and signal energies are accumulated separately
%   (numALS, denTrue) rather than stored as ratios. This keeps both the
%   classical Monte-Carlo NMSE, sum(num)/sum(den), and the per-trial
%   ratio, from which the paper's median and the NMSE>1 failure rates
%   are computed, available from one run.
%
% See also MAKE_PAPER_FIGURE, TRACE_DESIGN_TRAINING, TRACE_CHANNEL_ENERGIES.

clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
addpath(fileparts(thisDir)); setup_trace;

rng(31);

%% ---------------------------------------------------------------------
% Fixed system parameters
%% ---------------------------------------------------------------------
base.L   = 6;
base.Nr  = 8;
base.NRF = 4;
base.M   = 24;
base.P   = 8;
base.FEMbar_type = 'random-phase';
base.FRF_type    = 'random-phase';
base.min_angle_sep_deg = 5;

SNRdB   = 10;
Tvec    = [4 8 12 16 20 24];       % T = 4 is the zero-margin point
MC      = 200;
NSTARTS = 5;

thetaGrid = linspace(-pi/2,pi/2,1001);   % same grid for TRACE and SOMP

R  = base.L;
NT = numel(Tvec);

%% ---------------------------------------------------------------------
% Nested random-phase EM training family
%% ---------------------------------------------------------------------
parMaster = base;
parMaster.T = max(Tvec);

rng(1001);
[trainingMaster,diagMaster] = trace_design_training(parMaster);

trainingCell = cell(NT,1);
recCell      = cell(NT,1);
sompCell     = cell(NT,1);
muEM         = zeros(NT,1);

for it = 1:NT
    T = Tvec(it);

    tr = struct();

    % Nested EM codebook: same phases, more measurements as T increases.
    % The sqrt(T_max/T) rescaling keeps each EM column at unit expected
    % norm, so the SNR is comparable across T.
    tr.FEMbar = sqrt(parMaster.T/T) * trainingMaster.FEMbar(1:T,:);

    % Keep the RF combiner fixed for all T.
    tr.FRF = trainingMaster.FRF;

    trainingCell{it} = tr;
    recCell{it}  = trace_prepare_recovery(tr,base.Nr,thetaGrid);
    sompCell{it} = trace_somp_prepare(tr,base.Nr,thetaGrid);

    % Mutual coherence of eq. (16) for this codebook: the explanatory
    % variable for the T = 4 breakdown.
    F = tr.FEMbar ./ max(vecnorm(tr.FEMbar),eps);
    G = abs(F'*F);
    G(1:size(G,1)+1:end) = 0;
    muEM(it) = max(G(:));
end

fprintf('\nNMSE versus T at SNR = %.1f dB\n',SNRdB);
fprintf('Fixed RF minimum normalized gain = %.4f\n',diagMaster.min_RF_gain_rel);
for it = 1:NT
    fprintf('  T = %2d : mu(F_EMbar) = %.4f   (eq. 16)\n',Tvec(it),muEM(it));
end
fprintf('\n');

%% ---------------------------------------------------------------------
% Accumulated error and signal energies
%% ---------------------------------------------------------------------
numALS    = zeros(NT,MC);
numSOMP   = zeros(NT,MC);
numOracle = zeros(NT,MC);
denTrue   = zeros(NT,MC);
rankdefALS = false(NT,MC);

%% ---------------------------------------------------------------------
% Monte-Carlo loop
%% ---------------------------------------------------------------------
for mc = 1:MC

    % Generate ONE physical channel using the largest T. The same
    % H_EM, A(theta), alpha and pilots are reused for every T.
    parRef = base;
    parRef.T = parMaster.T;
    cleanMaster = trace_generate_clean_data(parRef,trainingMaster);

    for it = 1:NT

        T = Tvec(it);
        par = base;
        par.T = T;
        tr = trainingCell{it};

        % Rebuild the noiseless tensor for this T while preserving the
        % exact same physical channel realization: only T_H changes
        % (eq. (13)), because only F_EMbar changed.
        clean = cleanMaster;
        clean.par    = par;
        clean.FEMbar = tr.FEMbar;
        clean.FRF    = tr.FRF;
        clean.TH     = tr.FEMbar * clean.HEM;            % eq. (13)
        clean.TA     = tr.FRF * clean.A;                 % eq. (12)
        clean.Yclean = trace_cpdgen(clean.TA,clean.TH,clean.TS);   % eq. (11)

        data = trace_add_noise(clean,SNRdB);

        %% TRACE
        oa = struct();
        oa.init_method = 'svd';
        oa.nstarts     = NSTARTS;
        oa.maxit       = 200;
        oa.tol         = 1e-8;
        oa.verbose     = false;

        estALS  = trace_als(data.Y,R,oa);
        physALS = trace_physical_recovery(estALS,data,recCell{it});

        [numALS(it,mc),denTrue(it,mc)] = trace_channel_energies( ...
            data.A,data.HEM,data.alpha, ...
            physALS.A_hat,physALS.HEM_hat,physALS.alpha_hat);

        rankdefALS(it,mc) = physALS.rank_deficient_ZS;

        %% Joint SOMP
        physSOMP = trace_somp_baseline(data,sompCell{it},R);
        numSOMP(it,mc) = trace_channel_energies( ...
            data.A,data.HEM,data.alpha, ...
            physSOMP.A_hat,physSOMP.HEM_hat,physSOMP.alpha_hat);

        %% Oracle-support LS
        physOracle = trace_oracle_baseline(data);
        numOracle(it,mc) = trace_channel_energies( ...
            data.A,data.HEM,data.alpha, ...
            physOracle.A_hat,physOracle.HEM_hat,physOracle.alpha_hat);
    end

    if mod(mc,20)==0
        fprintf('  NMSE versus T: MC %d/%d\n',mc,MC);
    end
end

%% ---------------------------------------------------------------------
% Statistics
%% ---------------------------------------------------------------------
ratALS    = numALS    ./ max(denTrue,eps);   % per-trial NMSE
ratSOMP   = numSOMP   ./ max(denTrue,eps);
ratOracle = numOracle ./ max(denTrue,eps);

medALS    = median(ratALS,2);                % the paper reports the median
medSOMP   = median(ratSOMP,2);
medOracle = median(ratOracle,2);

% Classical Monte-Carlo NMSE (ratio of accumulated energies), kept for
% reference; the paper does not plot this one.
nmseALS    = sum(numALS,2)    ./ max(sum(denTrue,2),eps);
nmseSOMP   = sum(numSOMP,2)   ./ max(sum(denTrue,2),eps);
nmseOracle = sum(numOracle,2) ./ max(sum(denTrue,2),eps);

%% ---------------------------------------------------------------------
% Quick look
%% ---------------------------------------------------------------------
figure;
semilogy(Tvec,medALS,'-o','LineWidth',1.6,'MarkerSize',7); hold on;
semilogy(Tvec,medSOMP,'-.d','LineWidth',1.6,'MarkerSize',7);
semilogy(Tvec,medOracle,':^','LineWidth',1.8,'MarkerSize',7);

% Sufficient Kruskal threshold, eq. (29), for this configuration.
kA = min(base.NRF,R); kS = min(base.P,R);
Tmin = 2*R + 2 - kA - kS;
if exist('xline','file')
    xline(Tmin,':k','Kruskal bound','LabelVerticalAlignment','bottom');
else
    yl = ylim; plot([Tmin Tmin],yl,':k'); ylim(yl);
end

grid on; box on;
xlabel('Number of EM probing states, T'); ylabel('NMSE');
legend('TRACE-ALS','Joint SOMP','Oracle-support LS','Location','northeast');
title(sprintf('NMSE versus T at SNR = %g dB',SNRdB));
set(gca,'FontSize',11); set(gca,'XTick',Tvec);

%% ---------------------------------------------------------------------
% Save
%% ---------------------------------------------------------------------
outDir = fullfile(fileparts(thisDir),'results');
if ~exist(outDir,'dir'), mkdir(outDir); end
save(fullfile(outDir,'figure_c_data.mat'), ...
     'base','SNRdB','Tvec','MC','NSTARTS','R','muEM','Tmin', ...
     'numALS','numSOMP','numOracle','denTrue','rankdefALS', ...
     'medALS','medSOMP','medOracle','nmseALS','nmseSOMP','nmseOracle', ...
     'ratALS','ratSOMP','ratOracle');

fprintf('\nSaved %s\n',fullfile(outDir,'figure_c_data.mat'));

% Reproduce the numbers quoted in Section V for this panel.
t = @(v) find(Tvec==v,1);
fprintf('median NMSE at T=4        : %.3f\n',medALS(t(4)));
fprintf('NMSE>1 rate at T=4        : %.1f %%\n', ...
    100*mean(ratALS(t(4),:)>1));
fprintf('TRACE/oracle at T=4       : %.1fx\n',medALS(t(4))/medOracle(t(4)));
fprintf('T=4 -> T=8 improvement    : %.1fx\n',medALS(t(4))/medALS(t(8)));
fprintf('T=8 -> T=24 improvement   : %.1fx\n',medALS(t(8))/medALS(t(24)));
fprintf('EM coherence T=4 / T=24   : %.3f / %.3f\n',muEM(t(4)),muEM(t(24)));

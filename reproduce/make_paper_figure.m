function make_paper_figure(outfile)
%MAKE_PAPER_FIGURE Draw the three-panel results figure of the paper.
%
%   make_paper_figure
%   make_paper_figure('Fig_TRACE_results.pdf')
%
% Reads the three result files from results/ and renders the double-column
% figure used in the paper:
%
%   (a) channel NMSE versus SNR for three training dimension choices
%   (b) TRACE versus the greedy baseline and the oracle bound
%   (c) channel NMSE versus the number of EM probing states T
%
% All panels plot the MEDIAN over the Monte-Carlo realizations. The median
% is used rather than the mean because the distribution is heavy-tailed
% near the low-SNR threshold, where a handful of failed realizations
% (NMSE > 1) would otherwise dominate the average and hide the behaviour
% of the typical run. The per-realization matrices are kept in the saved
% files, so any other statistic can be computed without re-running.
%
% INPUT FILES (produced by the repro_fig_* scripts)
%   results/figure_a_data.mat   SNRdB, all_nmse{c}, labels
%   results/figure_b_data.mat   SNRdB, als, somp, oracle
%   results/figure_c_data.mat   Tvec, ratALS, ratSOMP, ratOracle, Tmin
%
% STYLE
%   Okabe-Ito colour-blind-safe palette with triple-redundant encoding:
%   every series differs in colour AND line style AND marker, so the
%   figure survives greyscale printing and the common forms of colour
%   vision deficiency.
%
% See also REPRO_FIG_A_TRAINING_DIMENSIONS, REPRO_FIG_B_BASELINES,
%          REPRO_FIG_C_PROBING_STATES.

if nargin < 1 || isempty(outfile)
    outfile = 'Fig_TRACE_results.pdf';
end

thisDir = fileparts(mfilename('fullpath'));
resDir  = fullfile(fileparts(thisDir),'results');

% Okabe-Ito: blue, vermillion, bluish green, dark yellow.
COLORS  = [0.000 0.447 0.698;
           0.835 0.369 0.000;
           0.000 0.620 0.451;
           0.702 0.478 0.000];
LINES   = {'-','--','-.',':'};
MARKERS = {'o','s','d','^'};
LW = 1.3; MS = 5; FS = 8;

fig = figure('Units','inches','Position',[1 1 7.16 2.45], ...
             'Color','w','PaperPositionMode','auto');

%% ------------------------------------------------------------ panel (a)
Fa = load(fullfile(resDir,'figure_a_data.mat'));
ax1 = subplot(1,3,1); hold(ax1,'on');
for c = 1:numel(Fa.all_nmse)
    y = median(Fa.all_nmse{c},2);
    plot(ax1,Fa.SNRdB,y,'Color',COLORS(c,:),'LineStyle',LINES{c}, ...
         'Marker',MARKERS{c},'LineWidth',LW,'MarkerSize',MS, ...
         'MarkerFaceColor','none');
end
finish_axes(ax1,'SNR (dB)','NMSE','(a) Training dimensions',FS);
legend(ax1,Fa.labels,'Location','southwest','FontSize',FS-1.5,'Box','off');

%% ------------------------------------------------------------ panel (b)
Fb = load(fullfile(resDir,'figure_b_data.mat'));
keep = Fb.SNRdB >= 1;                  % range shown in the paper
ax2 = subplot(1,3,2); hold(ax2,'on');
series = {median(Fb.als,2),median(Fb.somp,2),median(Fb.oracle,2)};
names  = {'TRACE','Joint SOMP','Oracle-support LS'};
for c = 1:3
    plot(ax2,Fb.SNRdB(keep),series{c}(keep), ...
         'Color',COLORS(c,:),'LineStyle',LINES{c},'Marker',MARKERS{c}, ...
         'LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor','none');
end
finish_axes(ax2,'SNR (dB)','NMSE','(b) Compressive configuration',FS);
legend(ax2,names,'Location','southwest','FontSize',FS-1.5,'Box','off');

%% ------------------------------------------------------------ panel (c)
Fc = load(fullfile(resDir,'figure_c_data.mat'));
ax3 = subplot(1,3,3); hold(ax3,'on');
seriesC = {median(Fc.ratALS,2),median(Fc.ratSOMP,2),median(Fc.ratOracle,2)};
for c = 1:3
    plot(ax3,Fc.Tvec,seriesC{c}, ...
         'Color',COLORS(c,:),'LineStyle',LINES{c},'Marker',MARKERS{c}, ...
         'LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor','none');
end
% The Kruskal bound of eq. (29): the T below which uniqueness is not
% guaranteed. The panel's point is that being just above it is not enough.
yl = ylim(ax3);
plot(ax3,[Fc.Tmin Fc.Tmin],yl,':k','LineWidth',0.9);
text(ax3,Fc.Tmin,yl(1)*2,' Kruskal bound','FontSize',FS-1.5, ...
     'VerticalAlignment','bottom');
ylim(ax3,yl);
finish_axes(ax3,'EM probing states, T','NMSE','(c) EM probing states',FS);
set(ax3,'XTick',Fc.Tvec);
legend(ax3,names,'Location','northeast','FontSize',FS-1.5,'Box','off');

%% ---------------------------------------------------------------- save
if exist('exportgraphics','file')==2
    exportgraphics(fig,outfile,'ContentType','vector');
else
    print(fig,outfile,'-dpdf','-painters');
end
fprintf('Wrote %s\n',outfile);

end

%% =====================================================================
function finish_axes(ax,xl,yl,ttl,FS)
set(ax,'YScale','log','FontSize',FS,'Box','on','XGrid','on','YGrid','on', ...
       'GridAlpha',0.15,'MinorGridAlpha',0.08,'Layer','top');
xlabel(ax,xl,'FontSize',FS);
ylabel(ax,yl,'FontSize',FS);
title(ax,ttl,'FontSize',FS,'FontWeight','normal');
end

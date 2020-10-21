clear variables
run('my_prefs.m')
loadpath = uigetdir('/Users/matthiasschickinger/PhD/TIRFM_Data/00_TIRFM_Analysis/NaCl_screens');
cd(loadpath)
load('lts_main_pop.mat')
close all
exppdf_mod = @(t,tau,Tmin,Tmax)1./tau.*exp(-t./tau)./(exp(-Tmin./tau)-exp(-Tmax./tau));
expcdf_mod = @(t,tau,Tmin,Tmax)(exp(-t./tau)-exp(-Tmin./tau))/(exp(-Tmax./tau)-exp(-Tmin./tau));
%colors = {[237 177 32]/255,[0 102 153]/255, [204 0 0]/255};
%cf = figure('Color', [1 1 1], 'Units', 'Pixels', 'Position', [500 500 128 115], 'PaperPositionMode', 'auto', ...
 %               'Visible','off');
[filename,savepath] = uiputfile('','',[paper_dir '/2017_Figures/NaClScreen/']);%inputdlg({'Enter filename (without extension):'},'Filename',1,{''});
cd(savepath)

%% Set limits, color, markers
BLIM = [10 50];
ULIM = [6.5 300];
% yellow-purple-green
%colors = {[.929 .694 .125],[.494 .184 .556],[.466 .674 .188]};
% azurblau-rot-orange
%colors = {[36 88 120]/255,[204 0 0]/255,[255 170 0]/255};
% blau verlauf
%colors = {[21 51 69]/255,[59 144 196]/255,[118 167 196]/255};
% lila verlauf
%colors = {[204 0 163]/255,[126 47 142]/255,[69 27 67]/255};
% yellow-orange verlauf 
%colors = {[247 127 21]/255,[.929 .694 .125],[247 226 21]/255};
%colors = {[230 94 43]/255,[255 170 0]/255,[.929 .694 .125]};
%colors = {[230 94 43]/255,[254 186 48]/255,[.929 .694 .125]};
% L2: cyan-aqua-darqua
%colors = {[128 255 255]./255,[0 128 255]./255,[0 0 128]./255}; mind = 1;
% L6: lemon-tangerine-maraschino
%colors = {[225 225 0]./255,[255 128 0]./255,[255 0 0]./255}; mind = 2;
% L10: magenta-plum
colors = {[255 0 255]./255,[128 0 128]./255}; mind = 3;
%cind = 2;
markers = {'^','o','s'};
%markers = {'x', 'd', '+'}
msize = [16,16,20];
% Scatter plot
close all
sz = [500 350];
sf = figure('Units','Points','Position',[1 1 sz],'PaperUnits','Points','PaperSize',sz,'PaperPositionMode','auto');
hold off
XYmin = [Inf Inf];
cind = 0;
for j = 1:numel(lifetimes)
    cind = cind+1;
    loglog(lifetimes{j}.MEAN(INdices{j},2),lifetimes{j}.MEAN(INdices{j},1),markers{mind},'Color',colors{cind},'MarkerSize',msize(mind),'LineWidth',.2,'MarkerFaceColor',[1 1 1]*1)
    hold on
    for k = 1:2
        XYmin(k) = min(XYmin(k),min(lifetimes{j}.MEAN(INdices{j},k)));
    end
end
XYmin = XYmin - [0.1 0.1];
ax = gca;
ax.Units = 'points';
ax.XLim = ULIM;
%xlim auto
ax.YLim = BLIM;
ax.TickDir = 'out';
%ax.FontName = 'Helvetica';
%ax.FontSize = 20;
%ax.YTick = [0 .5 1];
ax.XTick = [10,20,50,100,300];
ax.XMinorTick = 'on';
ax.XTickLabel = {};
% for i = 1:numel(ax.XTick)
%     ax.XTickLabel{i} = num2str(ax.XTick(i));
% end
ax.YTick = 10:10:80;%[10, 20, 40, 80];
ax.YMinorTick = 'off';
ax.YTickLabel = {};
% for i = 1:numel(ax.YTick)
%     ax.YTickLabel{i} = num2str(ax.YTick(i));
% end
ax.YMinorGrid = 'off';
%xlabel('lifetime (s)', 'FontSize', 14)
%ylabel('Cumulative frequency / Probability', 'FontSize', 14)
%l = legend({'data','MLE'},'Location','southeast', 'FontSize', 12);
box off
ax.Layer = 'bottom';
grid on
ax.GridLineStyle = '-';
ax.MinorGridLineStyle = '-';

%minimize whitespace (from Matlab help)
ax.Position(1) = 10;
ax.Position(2) = 10;
ax.Position(3:4) = sf.Position(3:4)-ax.Position(1:2)-[1 1];

%export_fig(filename(1:strfind(filename,'.')-1),'-eps','-r600', '-transparent')
print('-dpdf',[filename(1:strfind(filename,'.')-1) '.pdf'])

%% CDF plot unbound (state 2)
close all
sz = [100 100 320 200];
figure('Units','Points','Position',sz,'PaperUnits','Points','PaperSize',sz(3:4),'PaperPositionMode','auto')
hold off
cind = 0;
Tmax = 800;
tmax = 0;
for j = numel(lifetimes):-1:1
    cind = cind+1;
    lt = sort(vertcat(lifetimes{j}.ALL{INdices{j},2}));
    tmin = min(lt);
    %
    tmp = lt(floor(0.999*length(lt)));
    centers = logspace(-2,log10(Tmax),3e3);
    %centers = linspace(0,Tmax,3e3);
    cumcts = zeros(size(centers));
    for i = 1:length(centers)
        cumcts(i) = sum(lt<=centers(i));
    end
    cumcts = cumcts/cumcts(end);
    tmax = max(tmax,centers(find(cumcts>.999,1)));
    semilogx(centers,cumcts,'Color',colors{cind},'LineWidth',2)
    hold on
end
ax = gca;
ax.XLim = [tmin-0.1 tmax];
%xlim auto
ax.YLim = [0 1];
ax.YTick = [0 0.5 1];
ax.TickDir = 'out';
ax.FontName = 'Helvetica';
ax.FontSize = 12;
%ax.YTick = [0 .5 1];
ax.XTick = [.1, 1, 10, 100, 1000];
for i = 1:numel(ax.XTick)
    ax.XTickLabel{i} = num2str(ax.XTick(i));
end
%ax.YTickLabel = {};
%xlabel('lifetime (s)', 'FontSize', 14)
%ylabel('Cumulative frequency / Probability', 'FontSize', 14)
%l = legend({'data','MLE'},'Location','southeast', 'FontSize', 12);
box off
ax.Layer = 'top';
grid off

%minimize whitespace (from Matlab help)
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

%export_fig(filename(1:strfind(filename,'.')-1),'-eps','-r600', '-transparent')
print('-depsc','-r600','-tiff','-loose',[filename(1:strfind(filename,'.')-1) '_unbound.eps'])

%% CDF plot bound (state1)
close all
sz = [100 100 320 200];
figure('Units','Points','Position',sz,'PaperUnits','Points','PaperSize',sz(3:4),'PaperPositionMode','auto')
hold off
cind = 0;
Tmax = tmax;
for j = numel(lifetimes):-1:1
    cind = cind+1;
    lt = sort(vertcat(lifetimes{j}.ALL{INdices{j},1}));
    %Tmax = max(lt)+0.1;
    %
    tmp = lt(floor(0.999*length(lt)));
    centers = logspace(-2,log10(Tmax),3e3);
    %centers = linspace(0,Tmax,3e3);
    cumcts = zeros(size(centers));
    for i = 1:length(centers)
        cumcts(i) = sum(lt<=centers(i));
    end
    cumcts = cumcts/cumcts(end);
    %tmax = max(tmax,centers(find(cumcts>.999,1)));
    semilogx(centers,cumcts,'Color',colors{cind},'LineWidth',2)
    hold on
end
ax = gca;
ax.XLim = [tmin-0.1 tmax];
%xlim auto
ax.YLim = [0 1];
ax.YTick = [0 0.5 1];
ax.TickDir = 'out';
ax.FontName = 'Helvetica';
ax.FontSize = 12;
%ax.YTick = [0 .5 1];
ax.XTick = [.1, 1, 10, 100, 1000];
for i = 1:numel(ax.XTick)
    ax.XTickLabel{i} = num2str(ax.XTick(i));
end
%ax.YTickLabel = {};
%xlabel('lifetime (s)', 'FontSize', 14)
%ylabel('Cumulative frequency / Probability', 'FontSize', 14)
%l = legend({'data','MLE'},'Location','southeast', 'FontSize', 12);
box off
ax.Layer = 'top';
grid off

%minimize whitespace (from Matlab help)
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

%export_fig(filename(1:strfind(filename,'.')-1),'-eps','-r600', '-transparent')
print('-depsc','-r600','-tiff','-loose',[filename(1:strfind(filename,'.')-1) '_bound.eps'])

clear variables
run('my_prefs.m')
load('lts_main_pop.mat')
close all
exppdf_mod = @(t,tau,Tmin,Tmax)1./tau.*exp(-t./tau)./(exp(-Tmin./tau)-exp(-Tmax./tau));
expcdf_mod = @(t,tau,Tmin,Tmax)(exp(-t./tau)-exp(-Tmin./tau))/(exp(-Tmax./tau)-exp(-Tmin./tau));
colors = {[204 0 0]/255,[0 102 153]/255};
cf = figure('Color', [1 1 1], 'Units', 'Pixels', 'Position', [500 500 388 65], 'PaperPositionMode', 'auto', ...
                'Visible','off');
tmax = [1035 3825];
for j = 1:numel(lifetimes)
    for state = 1:2
        lt_state = sort(vertcat(lifetimes{j}.ALL{INdices{j},state}));
        Tmax = max(lt_state)+0.1;
        %
        tmp = lt_state(floor(0.999*length(lt_state)));
        centers = logspace(-2,log10(max(lt_state)),2e3);
        cumcts = zeros(size(centers));
        for i = 1:length(centers)
            cumcts(i) = sum(lt_state<=centers(i));
        end
        cumcts = cumcts/cumcts(end);
        %tmax = centers(find(cumcts>.999,1));
        %tmax = -tauhats(j,state)*log(0.01*(exp(-Tmin(state)/tauhats(j,state))-exp(-Tmax/tauhats(j,state))));

        % cdf
        subplot(1,2,state)
        hold off
        semilogx(centers,cumcts,'Color',[1 1 1]*.7,'LineWidth',1)
        hold on
        if centers(end)<tmax(state)
            semilogx([centers(end) tmax(state)],[1 1],'Color',[1 1 1]*.7,'LineWidth',1)
        end
        ts = logspace(log10(Tmin(state)),log10(max(centers(end),tmax(state))),2e3);
        semilogx(ts,expcdf_mod(ts,tauhats(j,state),Tmin(state),Tmax),'k--','LineWidth',.5)
        xlim([Tmin(state)-.1 max(Tmax,tmax(state))])
        ax = gca;
        ax.TickDir = 'out';
        ax.FontSize = 3;
        ax.YTick = [0 .5 1];
        for i = 1:numel(ax.XTick)
            ax.XTickLabel{i} = num2str(ax.XTick(i));
        end
        %ax.YTickLabel = {};
        %xlabel('lifetime (s)', 'FontSize', 14)
        %ylabel('Cumulative frequency / Probability', 'FontSize', 14)
        %l = legend({'data','MLE'},'Location','southeast', 'FontSize', 12);
        box off
        ax.Layer = 'top';
        %titlefontsize = 14;
        %title({['CDF for state ' num2str(state) ', \tau = ' num2str(round(tauhats(j,state),2)) ' s, N = ' num2str(Ns(j,state))],''}, 'FontSize', titlefontsize, 'Interpreter', 'tex')
    end
    %print('-dpng','-r150','Hist_KDE_CDF.png')
    %print('-depsc','-r0','-tiff','Hist_KDE_cdf.eps')
end
[filename,pathname] = uiputfile;%inputdlg({'Enter filename (without extension):'},'Filename',1,{''});
cd(pathname)
export_fig(filename(1:strfind(filename,'.')-1),'-eps','-png','-r600', '-transparent')
%export_fig([filename '.png'], '-r600')
clear variables
run('my_prefs.m')
loadpath = uigetdir('/Users/matthiasschickinger/PhD/TIRFM_Data/00_TIRFM_Analysis/2017_07_06_var8_mainpops');
cd(loadpath)
load('lts_main_pop.mat')
close all
exppdf_mod = @(t,tau,Tmin,Tmax)1./tau.*exp(-t./tau)./(exp(-Tmin./tau)-exp(-Tmax./tau));
expcdf_mod = @(t,tau,Tmin,Tmax)(exp(-t./tau)-exp(-Tmin./tau))/(exp(-Tmax./tau)-exp(-Tmin./tau));
colors = {[204 0 0]/255,[0 102 153]/255};
cf = figure('Color', [1 1 1], 'Units', 'Pixels', 'Position', [500 500 128 115], 'PaperPositionMode', 'auto', ...
                'Visible','off');
[filename,savepath] = uiputfile('','',[paper_dir '/2017_Figures/8merVariants/']);%inputdlg({'Enter filename (without extension):'},'Filename',1,{''});
cd(savepath)
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
        tmax = centers(find(cumcts>.999,1));

        % cdf
        hold off
        semilogx(centers,cumcts,'Color',[1 1 1]*.7,'LineWidth',2)
        hold on
        ts = logspace(log10(Tmin(state)),log10(centers(end)),2e3);
        semilogx(ts,expcdf_mod(ts,tauhats(j,state),Tmin(state),Tmax),'k--','LineWidth',.75)
        xlim([Tmin(state)-.1 tmax])
        ax = gca;
        ax.TickDir = 'out';
        ax.FontName = 'Helvetica Neue';
        ax.FontSize = 12;
        ax.YTick = [0 .5 1];
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
        %titlefontsize = 14;
        %title({['CDF for state ' num2str(state) ', \tau = ' num2str(round(tauhats(j,state),2)) ' s, N = ' num2str(Ns(j,state))],''}, 'FontSize', titlefontsize, 'Interpreter', 'tex')
        export_fig([filename(1:strfind(filename,'.')-1) 'state' num2str(state)],'-eps','-r600', '-transparent')
    end
    %print('-dpng','-r150','Hist_KDE_CDF.png')
    %print('-depsc','-r0','-tiff','Hist_KDE_cdf.eps')
end
%export_fig([filename '.png'], '-r600')
for i = 1:length(paths)
cd(paths{i})
path_out = cd;
mkdir([path_out filesep 'positions'])

clear data
%load('all_data_struct.mat')
load('data_spot_pairs.mat', 'data')

close all
ps = figure('Position', [0 0 800 800], 'Visible', 'Off');
colors = {'r';'g'};
for m = 1:size(data,1)
    for s = 1:size(data{m},1)
        for ch = 1:2
            plot_pos0 = data{m}{s,ch}.pos0(data{m}{s,ch}.vwcm.pos(:,1)>0,:);
            plot_pos = data{m}{s,ch}.vwcm.pos(data{m}{s,ch}.vwcm.pos(:,1)>0,:);
            plot_means = data{m}{s,ch}.vwcm.means100;
            subplot(2,2,ch)
            hold off
            plot(plot_pos(:,1),plot_pos(:,2),[colors{ch} '.'],'MarkerSize',8)
            hold on
            plot(plot_means(:,1),plot_means(:,2),'Color', [1 1 1]*.7,'LineWidth', 2)
            plot(plot_pos0(:,1),plot_pos0(:,2),'k-x','MarkerSize', 16, 'LineWidth', 2)
            
            axis equal
            if size(plot_pos0,1)>0
            xlim([min(plot_pos0(:,1))-5 max(plot_pos0(:,1))+5])
            ylim([min(plot_pos0(:,2))-5 max(plot_pos0(:,2))+5])
            end
            grid on
            set(gca, 'YDir', 'normal')
            title(['Pos and Pos0, channel ' num2str(ch)])

            subplot(2,2,ch+2)
            hold off
            plot(plot_pos(:,1)-plot_pos0(:,1),plot_pos(:,2)-plot_pos0(:,2),[colors{ch} '.'],'MarkerSize',8)
            axis equal
            xlim([-5 5])
            ylim([-5 5])
            grid on
            set(gca, 'YDir', 'normal', 'XTick', -5:1:5, 'YTick', -5:1:5)            
            title(['Deviation from Pos0, channel ' num2str(ch)])
        end
    % Save as image
    print('-dpng', '-r96', [path_out filesep 'positions' filesep 'positions_m' num2str(m) '_s' num2str(s) '.png'])
    end
end
end
display('Done')
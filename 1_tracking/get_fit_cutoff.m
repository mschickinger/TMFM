function [ fit_cutoff ] = get_fit_cutoff( mov, cut, merged_itraces, ...
                            avg_img, channel, N_movie, m, r_integrate )
% Determine intensity or frame cutoff value for pos_in_frame assignment for
% one movie in one channel. Function called by TWO_COLOR_TRACKER from v13.

scrsz = get(0,'ScreenSize');
figure('Position', scrsz);

% sort spots by iEndval
iEndval = zeros(size(merged_itraces,1),1);
for j=1:length(merged_itraces)
    iEndval(j) = median(merged_itraces{j}(end-500:end,4));
end   
[tmp_val, tmp_spot] = sort(iEndval);
iEndval_sorted = [tmp_spot tmp_val];
levels = [mean(iEndval) median(iEndval_sorted(1:ceil(length(iEndval_sorted)/10),2))];

%create fit_cutoff array
if cut == 1
    def_fc = iEndval_sorted(1,2);
    fit_cutoff = zeros(size(merged_itraces,1),1);
elseif cut == 2
    def_fc = length(mov.frames);
    fit_cutoff = def_fc.*ones(size(merged_itraces,1),1);
end

%cycle through spots, in iEndval ascending order
for j = 1:size(merged_itraces,1)
    act_spotnum = iEndval_sorted(j,1);
    x_0 = round(merged_itraces{act_spotnum}(1,2));
    y_0 = round(merged_itraces{act_spotnum}(1,3));
    x_0_end = round(merged_itraces{act_spotnum}(end,2));
    y_0_end = round(merged_itraces{act_spotnum}(end,3));

    subplot('Position', [0.05,0.55,0.65,0.4])
    hold off
    plot (iEndval_sorted(:,2), [channel '.'], 'MarkerSize', 8);
    hold on
    plot (j,iEndval_sorted(j,2), 'ko', 'MarkerSize', 10);
    for l = 1:length(levels)
        plot ([1 length(iEndval_sorted)],[1 1]*levels(l),'-');
    end
    title('Average end value distribution (ascending order)')

    subplot('Position', [0.75,0.55,0.25,0.4])
    hold off
    plot_subframe(avg_img{1}, x_0, y_0, 6)
    hold on
    ellipse(r_integrate, r_integrate, 0, x_0, y_0, channel)
    title('Averaged over first 100 frames');
    axis square
    set(gca, 'YDir', 'normal')

    subplot('Position', [0.75,0.05,0.25,0.4])
    hold off
    plot_subframe(avg_img{2}, x_0_end, y_0_end, 6)
    hold on
    ellipse(r_integrate, r_integrate, 0, x_0_end, y_0_end, channel)
    title('Averaged over last 100 frames');
    axis square
    set(gca, 'YDir', 'normal')

    subplot('Position', [0.05,0.05,0.65,0.4])
    hold off
    plot(merged_itraces{act_spotnum}(:,1), ...
        merged_itraces{act_spotnum}(:,4),...
        ['-' channel(1)], 'LineWidth', 0.5)
    hold on
    plot(merged_itraces{act_spotnum}(:,1), ...
        merged_itraces{act_spotnum}(:,5),...
        '-k', 'LineWidth', 0.25)
    tmp = ylim;
    if cut == 1
        plot(merged_itraces{act_spotnum}(:,1), ones(1,size(merged_itraces{act_spotnum},1)).*def_fc, ...
            'color', [1 1 1].*0.7, 'LineStyle', '-.')
    elseif cut == 2
        plot([def_fc def_fc], [tmp(1) tmp(2)],'color', [1 1 1].*0.7, 'LineStyle', '-.')
    end
    ylim(tmp)
    set(gca, 'YTick', tmp(1):500:tmp(2), 'Layer', 'top')
    grid on
    xlim([merged_itraces{act_spotnum}(1,1) merged_itraces{act_spotnum}(end,1)])
    title(['Intensity trace: cutoff selection for spots in movie #' num2str(m) ' of ' num2str(N_movie)])

    % determine cutoff
    if cut == 1
        tmp = inputdlg(['Cutoff intensity for movie #' num2str(m) ', ' channel ' channel, spot #' ...
            num2str(act_spotnum)], 'Cutoff intensity', 1, {num2str(def_fc)});
        if isempty(tmp)
            break
        end
        fit_cutoff(act_spotnum) = str2double(tmp);
        def_fc = (str2double(tmp)>0)*str2double(tmp);
    elseif cut == 2
        h = impoint(gca);
        if size(h,1) == 0
            break
        end
        tmp = getPosition(h);
        tmp = tmp(1);
        fit_cutoff(act_spotnum) = tmp;
        def_fc = (tmp>0)*tmp;
    end        
end

close all

end


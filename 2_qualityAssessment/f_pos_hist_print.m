function f_pos_hist_print(path_out, varargin)

%% parse input

p = inputParser;

addRequired(p, 'path_out', @isdir);
addParameter(p, 'filter', 'RMS', @ischar);
addParameter(p, 'method', 'vwcm', @ischar);
addParameter(p, 'YLIM', [0 3], @isvector);
addParameter(p, 'horzAx', 'frames', @ischar);

parse(p, path_out, varargin{:});
YLIM = p.Results.YLIM;

%% load data
clear avg_img r_integrate channel chb cut fit_cutoff data

cd(path_out)
mkdir([path_out filesep p.Results.method '_traces']);

load('data_proc.mat');
load('data_archive.mat', 'avg_img', 'r_integrate');
load('data_plot.mat');
load('data_spot_pairs.mat', 'data');

disp('All data loaded.')

%% Set parameters
N_movie = size(data,1);

if chb == 1
    chm = 2;
else
    chm = 1;
end

type = cell(2,1);
type{chm} = 'mobile';
type{chb} = 'fixed';

if ~exist('cut', 'var')
    cut = questdlg('Intensity threshold or maximum frame?','Cutoff method','Intensity','Frame','Frame');
end
if ischar(cut)
    if strcmp(cut, 'Intensity')
        cut_plot = [0 1];
    elseif strcmp(cut,'Frame')
        cut_plot = [1 0];
    end
else
    if cut == 1
        cut_plot = [0 1];
    elseif cut == 2
        cut_plot = [1 0];
    end
end

if strcmp(p.Results.horzAx, 'time')
    if ~exist('time_per_frame', 'var')
        prompt = cell(N_movie, 1);
        default_ans = prompt;
        for m = 1:N_movie
            prompt{m} = ['Enter time per frame in movie #' num2str(m) ':'];
            default_ans{m} = '100';
        end
        tmp = inputdlg(prompt, 'Time per frame', 1, default_ans);
        time_per_frame = zeros(N_movie,1);
        for m = 1:N_movie
            time_per_frame(m) = str2double(tmp{m});
        end
    end
end

%% Plot and print
for m=1:N_movie
    switch p.Results.horzAx
        case 'frames'
            if ceil(size(data{m}{1,1}.vwcm.pos,1)/1000) == 1
                XTickDiv = 20;
            elseif ceil(size(data{m}{1,1}.vwcm.pos,1)/1000) <= 10
                XTickDiv = 100;
            else
                XTickDiv = 1000;
            end
        case 'time'
            tpf_factor = time_per_frame(m)/500;
            if ceil(size(data{m}{1,1}.vwcm.pos,1)*tpf_factor/120) == 1
                XTickDiv = 12/tpf_factor;
                timestring = 'seconds';
            else
                if ceil(size(data{m}{1,1}.vwcm.pos,1)*tpf_factor/120) <= 10
                    XTickDiv = 60/tpf_factor;
                else
                    XTickDiv = 240/tpf_factor;
                end
                timestring = 'minutes';
            end
    end
    XTC = cell(1,floor(size(data{m}{1,1}.vwcm.pos,1)/XTickDiv));
    switch p.Results.horzAx
        case 'frames'
            for i = 1:floor(size(data{m}{1,1}.vwcm.pos,1)/(5*XTickDiv))
                XTC{5*i+1} = i*5*XTickDiv;
            end
        case 'time'
            for i = 1:floor(size(data{m}{1,1}.vwcm.pos,1)/(5*XTickDiv))
                XTC{5*i+1} = i*5*XTickDiv*tpf_factor/60;
            end
    end
    for s=1:size(data{m},1)
        figure('Visible','off', 'PaperPositionMode', 'manual','PaperUnits',...
        'centimeters', 'PaperPosition', [0 0 29.7 21]);
        for c=1:2
            switch p.Results.method
                
                case 'vwcm'
                    plot_data = data{m}{s,c}.vwcm; %adjust method for switching between FITTED and VWCM traces!!!
                case 'gF'
                    plot_data = data{m}{s,c}.gF;
            end
            tmp_plot = (1:size(plot_data.pos,1))';
            tmp_plot = [tmp_plot(plot_data.pos(:,1)>0) plot_data.r(plot_data.pos(:,1)>0)];
            %tmp_rms = [plot_data.frames(plot_data.frames>0) plot_data.rms10(plot_data.frames>0)];
                
            % Plot some subplots (first row for channel c)
                subplot('Position', [0.045 ((2-c)*0.5+0.3) 0.65 0.165])
                hold off
                if ~isempty(tmp_plot(tmp_plot(:,2)>5,1))
                    for i = tmp_plot(tmp_plot(:,2)>5,1)
                        plot([i i], YLIM, '-', 'Color', [1 1 1].*.7);
                    hold on
                    end
                end
                plot(tmp_plot(:,1), tmp_plot(:,2),['.' channel{c}(1)], 'MarkerSize', 5);
                hold on
                switch p.Results.filter
                    case 'RMS'
                        plot(tmp_plot(:,1),plot_data.rms10(plot_data.pos(:,1)>0),'-k', 'LineWidth', 0.1);
                    case 'median'
                        plot(tmp_plot(:,1),medfilt1_trunc(plot_data.r(plot_data.pos(:,1)>0),11),'-c', 'LineWidth', 0.1);
                end
                xlim([0 size(plot_data.pos,1)]) %DO PROPERLY IN THE LONG RUN!!!!
                ylim(YLIM) %WATCH OUT - NEED TO ADJUST FOR OTHER PLOTS !!!!!
                set(gca, 'XTick', 0:XTickDiv:size(plot_data.pos,1), 'XTickLabel', XTC, 'TickDir', 'in', 'Layer', 'top');
                grid on 
                title([channel{c} ' (' type{c} ') spot: Excursion / Rms10 traces. Spot pair # ' num2str(s) ' of '...
                    num2str(size(data{m},1)) ' / movie: ' num2str(m) ' of ' num2str(N_movie)])
                ylabel('px')
                switch p.Results.horzAx
                    case 'frames'
                        xlabel('frames')
                    case 'time'
                        xlabel(['time (' timestring ')'])
                end

                
                subplot('Position', [0.70 ((2-c)*0.5+0.3) 0.0625 0.165])
                hold off
                if sum(plot_data.pos(:,1) > 0) > 0
                    hist(plot_data.r(plot_data.r<=5),0:.05:5);
                else
                    hist(0,0:.05:5);
                end
                h1 = findobj(gca, 'Type', 'patch');
                set(h1, 'FaceColor', channel{c}(1), 'EdgeColor', channel{c}(1));
                hold on
                set(gca,'view',[90 -90], 'XTickLabel', {},'YTickLabel', {}, 'TickDir', 'in', 'Layer', 'top');
                %title('Histogram excursion')
                xlim(YLIM)

                
                subplot('Position', [0.7675 ((2-c)*0.5+0.3) 0.0625 0.165])
                hold off
                if sum(plot_data.pos(:,1) > 0) > 0
                    hist(plot_data.rms10(plot_data.rms10<=5),0:.05:5);
                else
                    hist(0,0:.05:5);
                end
                h2 = findobj(gca, 'Type', 'patch');
                set(h2, 'FaceColor', 'k');
                hold on
                set(gca,'view',[90 -90], 'XTickLabel', {},'YTickLabel', {},'XAxisLocation', 'bottom', 'TickDir', 'in', 'Layer', 'top');
                %title('Histogram RMS10 deviation')
                xlim(YLIM)

                
                subplot('Position', [0.83 ((2-c)*0.5+0.3) 0.165 0.165])
                hold off
                plot(0, 0, [channel{c}(1) 'x'])
                plot(plot_data.dispmed101(plot_data.r>5,1),plot_data.dispmed101(plot_data.r>5,2), 'o', 'Color', [1 1 1]*.7, 'MarkerSize', 3);
                hold on
                plot(plot_data.dispmed101(plot_data.r<=5,1),plot_data.dispmed101(plot_data.r<=5,2), [channel{c}(1) '+'], 'MarkerSize', 3);
                xlim([-5 5])
                ylim([-5 5])
                TC = cell(1,11); TC{1} = -5; TC{6} = 0; TC{11} = 5;
                set(gca, 'XTick', -5:5, 'YTick', -5:5, 'XTickLabel', TC, 'YTickLabel', TC);
                grid on, axis square
                title('Positions (drift corrected)')

                
                % Plot more subplots (second row for channel c)
                subplot('Position',[0.045, ((2-c)*0.5+0.05), 0.65, 0.175])
                hold off
                plot(data{m}{s,c}.itrace, ['-' channel{c}(1)], 'LineWidth', 0.75);
                hold on
                plot(data{m}{s,c}.med_itrace, 'k-', 'Linewidth', 0.5);
                plot(cut_plot*[fit_cutoff{m,c}(s) fit_cutoff{m,c}(s) ; 0 length(data{m}{s,c}.itrace)],...
                cut_plot*[0 max(data{m}{s,c}.itrace) ; fit_cutoff{m,c}(s) fit_cutoff{m,c}(s)],'k-');
                xlim([0 size(data{m}{s,c}.itrace,1)])
                ylim([.9*min(data{m}{s,c}.itrace) max([data{m}{s,c}.itrace; 1])])
                set(gca, 'XTick', 0:XTickDiv:size(data{m}{s,c}.itrace,1), 'XTickLabel', XTC, 'TickDir', 'in', 'Layer', 'top');
                grid on
                switch p.Results.horzAx
                    case 'frames'
                        title(['Summed pixel Intensities. Cutoff: ' num2str(double(fit_cutoff{m,c}(s)))])
                        xlabel('frames')
                    case 'time'
                        title(['Summed pixel Intensities. Cutoff: ' num2str(double(fit_cutoff{m,c}(s)*tpf_factor/60)) ' ' timestring '.'])
                        xlabel(['time (' timestring ')'])
                end

                
                subplot('Position',[0.70, ((2-c)*0.5+0.04), 0.155, 0.175])
                hold off
                x_0 = round(mean(data{m}{s,c}.pos0(1:100,1)));
                y_0 = round(mean(data{m}{s,c}.pos0(1:100,2)));
                if size(avg_img,2) == 4
                    plot_subframe(avg_img{m,c}, x_0, y_0, 12);
                elseif size(avg_img,2) == 2
                    plot_subframe(avg_img{m,c}{1}, x_0, y_0, 12);
                end
                hold on        
                ellipse(r_integrate, r_integrate, 0, x_0, y_0, channel{c});
                title(['Start Pos' sprintf('\n') num2str(x_0) ', ' num2str(y_0)]);
                axis square
                set(gca, 'YDir', 'normal', 'XTickLabel', {}, 'YTickLabel', {});

                
                subplot('Position',[0.85, ((2-c)*0.5+0.04), 0.155, 0.175])
                hold off
                x_0 = round(mean(data{m}{s,c}.pos0(end-99:end,1)));
                y_0 = round(mean(data{m}{s,c}.pos0(end-99:end,2)));
                if size(avg_img,2) == 4
                    plot_subframe(avg_img{m,c+2}, x_0, y_0, 12);
                elseif size(avg_img,2) == 2
                    plot_subframe(avg_img{m,c}{end}, x_0, y_0, 12);
                end
                hold on        
                ellipse(r_integrate, r_integrate, 0, x_0, y_0, channel{c});
                title(['End Pos' sprintf('\n') num2str(x_0) ', ' num2str(y_0)]);
                axis square
                set(gca, 'YDir', 'normal', 'XTickLabel', {}, 'YTickLabel', {});
        end
        %%
        % Save as image
        %set(gcf, 'PaperPositionMode', 'auto')
        print('-dpng', '-r96', [path_out filesep p.Results.method '_traces' filesep 'traces_RMS_hist_m' num2str(m) '_s' num2str(s) '.png']);
        close(gcf)
    end
    display(['Done printing figures for movie #' num2str(m)]) 
end
end


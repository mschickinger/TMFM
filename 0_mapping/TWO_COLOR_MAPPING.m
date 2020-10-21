%% include TRE control and nn optimisation

%% STARTUP
clc, clear all, close all
path0 = cd; 
run('my_prefs')
%% choose colors
channel = cell(2,1);
channel{1} = 'red';
channel{2} = 'green';

%% LOAD STACK OF MOVIES/FILES
pname=uigetdir(mapping_dir,'Choose the folder with all .fits files for BEADS mapping');
files_ch1 = pickFirstFitsFiles(pname, channel{1});
files_ch2 = pickFirstFitsFiles(pname, channel{2});

N_movie = length(files_ch1); % number of movies
if length(files_ch1) ~= length(files_ch2)
    disp('WARNING: not same number of movie files!')
end

path_out = [pname filesep datestr(now, 'yyyy-mm-dd_HH-MM') '_BEADS mapping'];
mkdir(path_out)

%% SET PARAMETER
input = {'First Frame:', 'Last Frame (-1=all):', ['Sequence ' channel{1} ':'], ['Sequence ' channel{2} ':'],... % sample options
    'Radius of peak [pixel]:', 'Number of frames for average [frames]:'};
input_default = {'1', '-1', '1', '1', '4', '1'};
tmp = inputdlg(input, 'Parameters', 1, input_default);
first = round(str2double(tmp(1))); % first image to read from file
last = round(str2double(tmp(2))); % last image to read from file
%determine sequences 
sequence_ch1 = zeros(1, size(tmp{3},2));
for i=1:size(tmp{3},2)
    if(tmp{3}(i) == '1')
        sequence_ch1(1,i) =1;
    end
end
sequence_ch2 = zeros(1, size(tmp{4},2));
for i=1:size(tmp{4},2)
    if(tmp{4}(i) == '1')
        sequence_ch2(1,i) =1;
    end
end
r_find = str2double(tmp(5)); % radius used to find spots and spots pairs
N_frames_for_average = str2double(tmp(6)); 



%% generate movie classes
ch1 = cell(N_movie,1);
ch2 = cell(N_movie,1);  

for i=1:N_movie
    ch1{i} = movie(pname, files_ch1{i}, first, last, sequence_ch1); % pname, fname, first, last, sequence
    ch2{i} = movie(pname, files_ch2{i}, first, last, sequence_ch2); % pname, fname, first, last, sequence
end


%% determine thresholds and find peaks
peaks_raw = zeros(0,5);

for i=1:N_movie
    [~,p1] = ch1{i}.get_h_min(r_find, N_frames_for_average, 'autosigma', 2);
    [~,p2] = ch2{i}.get_h_min(r_find, N_frames_for_average, 'autosigma', 2);
    
    % map peaks
    trace_map = map_traces(p1(:,1:2), p2(:,1:2), p2(:,1:2), r_find*2)+1; %map the traces from average positions

    tmp = zeros(size(trace_map,1),5);
    % combine pairs
    for j=1:size(trace_map,1)
        tmp(j,:) = [p1(trace_map(j,1), 1:2)+1 p2(trace_map(j,2), 1:2)+1 i]; %x_1 y_1 x_2 y_2 frame
    end
    
    peaks_raw = [peaks_raw; tmp];
end

N_peaks_raw = size(peaks_raw,1);
display(['You have ' num2str(N_peaks_raw) ' pairs'])

close all

%% compute averages images
avg_img = cell(N_movie, 2);

for i=1:N_movie
    avg_img{i, 1} = ch1{i}.average_image(N_frames_for_average);
    avg_img{i, 2} = ch2{i}.average_image(N_frames_for_average);
end

%% correction of average images and peak positions in case of automapping
if strcmp(questdlg('Images generated with automapping?', 'Automapping?', 'Yes'), 'Yes')
    for i = 1:N_movie
        for ch = 1:2
            avg_img{i,ch} = rot90(avg_img{i,ch},3);
        end
    end
    tmp = [(-peaks_raw(:,2)+ch1{1}.sizeY+1) peaks_raw(:,1) (-peaks_raw(:,4)+ch1{1}.sizeY+1) peaks_raw(:,3)];
    peaks_raw(:,1:4) = tmp;
    display('Rotated average images and peak positions 90 degrees.')
end

%% Fit psf to spots_ s_x ~ s_y

sigma_init = 1.5;
w_fit= 4;

ch1_fit_raw = zeros(N_peaks_raw, 7); 
ch1_fit_err_raw = zeros(N_peaks_raw, 7); 
ch2_fit_raw = zeros(N_peaks_raw, 7); 
ch2_fit_err_raw = zeros(N_peaks_raw, 7); 

h = waitbar(0,'Fitting spots... please wait');

for i=1:N_peaks_raw
    %display(['Fitting spot ' num2str(i) ' of ' num2str(N_peaks_raw)])
    
    % channel 1
    x1 = round(peaks_raw(i,1));
    y1 = round(peaks_raw(i,2));
    [c, c_err, ci, area] = fit_gauss2d_mainaxis_bg(x1, y1,sigma_init , w_fit, avg_img{peaks_raw(i, 5),1});
    ch1_fit_raw(i,:) = c;
    ch1_fit_err_raw(i,:) = c_err;
    
    % channel 2
    x2 = round(peaks_raw(i,3));
    y2 = round(peaks_raw(i,4));
    [c, c_err, ci, area] = fit_gauss2d_mainaxis_bg(x2, y2, sigma_init , w_fit, avg_img{peaks_raw(i, 5),2});
    ch2_fit_raw(i,:) = c;
    ch2_fit_err_raw(i,:) = c_err;
    
    waitbar( i/N_peaks_raw , h, ['Fitting spot... ' num2str(i) ' of ' num2str(N_peaks_raw) ' done']) % update waitbar

end

close(h)


%% SORT OUT #1

criteria = ones(N_peaks_raw,4 );

criteria(:,1:2) = filter_spots(ch1_fit_raw(:,3:4), [0.8 1.25], [1 1]);
criteria(:,3:4) = filter_spots(ch2_fit_raw(:,3:4), [0.8 1.25], [1 1]);
accepted = [criteria(:,1) & criteria(:,2) & criteria(:,3) & criteria(:,4)];

%remove not-accepted spots
ch1_fit = ch1_fit_raw(accepted==1, :);
ch1_fit_err = ch1_fit_err_raw(accepted==1, :);
ch2_fit = ch2_fit_raw(accepted==1, :);
ch2_fit_err = ch2_fit_err_raw(accepted==1, :);
peaks_raw2 = peaks_raw(accepted==1, :);
peaks_raw2 = [ch1_fit(:,1:2) ch2_fit(:,1:2) peaks_raw2(:,5)]; % use fitted poistions for further analysis

%peaks = peaks_raw(accepted==1, :);

%peaks = [ch1_fit(:,1:2) ch2_fit(:,1:2) peaks(:,5)]; % use fitted poistions for further analysis

%%
plot_discarded = 0; %strcmp(questdlg('Plot discarded spots?','Plot discarded','Yes','No','No'), 'Yes');
if plot_discarded
    path_out_discarded = [path_out filesep 'discarded_selection_n1'];
    mkdir(path_out_discarded)
end

close all
fig_dim =1*[20 10];
cur_fig = figure('Visible','off', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
colormap gray
w_plot = 10;

%
if plot_discarded
    display('Discarding spots...')
    for i=1:N_peaks_raw
        if ~accepted(i) % discarded spot
            message = {'Sigma ratio red: OK', 'Sigma ratio green: OK', 'Spotsize red: OK', 'Spotsize green: OK'};
            if criteria(i,1)==0
                message{1} = 'Sigma ratio red: BAD';
            end
            if criteria(i,2)==0
                message{2} = 'Sigma ratio green: BAD';
            end
            if criteria(i,3)==0
                message{3} = 'Spotsize red: BAD';
            end
            if criteria(i,4)==0
                message{4} = 'Spotsize green: BAD';
            end
            x_1 = ch1_fit_raw(i,1);
            y_1 = ch1_fit_raw(i,2);

            x_2 = ch2_fit_raw(i,1);
            y_2 = ch2_fit_raw(i,2);

            subplot(1, 2, 1)
            plot_subframe(avg_img{peaks_raw(i, 5), 1}, x_1, y_1, w_plot), hold on
            plot(x_1, y_1, 'g.')
            ellipse(ch1_fit_raw(i,3), ch1_fit_raw(i,4), -ch1_fit_raw(i,5), x_1, y_1, 'g'  );
            title({['Pair ' num2str(i) ' of '  num2str(N_peaks_raw) ' at (' num2str(round(x_1)) ',' num2str(round(y_1)) ') in red channel'], message{1}, message{3}})
            axis square
            hold off

            subplot(1, 2, 2)
            plot_subframe(avg_img{peaks_raw(i, 5), 2}, x_1, y_1, w_plot), hold on
            plot(x_2, y_2, 'r.')
            ellipse(ch2_fit_raw(i,3), ch2_fit_raw(i,4), -ch2_fit_raw(i,5), x_2, y_2, 'r'  );
            title({['Pair ' num2str(i) ' of '  num2str(N_peaks_raw) ' at (' num2str(round(x_2)) ',' num2str(round(y_2)) ') in green channel'], message{2}, message{4}})
            axis square
            hold off

            print(cur_fig, '-dtiff', '-r150',  [path_out_discarded filesep 'Discarded_' num2str(i) '.tif'])
        end 
    end
end

display(['Discarded ' num2str(sum(~accepted)) ' spots.'])
close all

N_pairs = size(peaks_raw2,1);


%% Fit#2: psf to spots s_x = s_y

sigma_init = 1.5;
w_fit= 4;

ch1_2fit= zeros(N_pairs, 5); 
ch1_2fit_err = zeros(N_pairs, 5); 
ch2_2fit = zeros(N_pairs, 5); 
ch2_2fit_err = zeros(N_pairs, 5); 

h = waitbar(0,'Fitting spots... please wait');

for i=1:N_pairs
    %display(['Fitting spot ' num2str(i) ' of ' num2str(N_pairs)])
    
    % channel 1
    x1 = round(ch1_fit(i,1));
    y1 = round(ch1_fit(i,2));
    [c, c_err, ci, area] = fit_gauss2d_bg(x1, y1,sigma_init , w_fit, avg_img{peaks_raw2(i, 5),1});

    ch1_2fit(i,:) = c;
    ch1_2fit_err(i,:) = c_err;
    
  
    %channel 2
    x2 = round(ch2_fit(i,1));
    y2 = round(ch2_fit(i,2));
    [c, c_err, ci, area] = fit_gauss2d_bg(x2, y2, sigma_init , w_fit, avg_img{peaks_raw2(i, 5),2});
    ch2_2fit(i,:) = c;
    ch2_2fit_err(i,:) = c_err;
    
    waitbar( i/N_pairs , h, ['Fitting spot... ' num2str(i) ' of ' num2str(N_pairs) ' done']) % update waitbar

end

close(h)


peaks = [ch1_2fit(:,1:2) ch2_2fit(:,1:2) peaks_raw2(:,5)]; % use fitted poistions for further analysis

N_peaks= size(peaks,1);
%% SORT OUT: remove spots which are too large
%{
criteria = ones(N_peaks_raw,2 );
criteria(:,1) = filter_spots_v02(ch1_fit(:,3), [1 1.5]);
criteria(:,2) = filter_spots_v02(ch2_fit(:,3), [1 1.5]);
accepted = [criteria(:,1) & criteria(:,2)];
%%
%remove not-accepted spots
ch1_fit = ch1_fit_raw(accepted==1, :);
ch1_fit_err = ch1_fit_err_raw(accepted==1, :);
ch2_fit = ch2_fit_raw(accepted==1, :);
ch2_fit_err = ch2_fit_err_raw(accepted==1, :);
peaks = peaks_raw(accepted==1, :);

peaks = [ch1_fit(:,1:2) ch2_fit(:,1:2) peaks(:,5)]; % use fitted poistions for further analysis


%%
plot_discarded = strcmp(questdlg('Plot discarded spots?','Plot discarded','Yes','No','No'), 'Yes');
if plot_discarded
    path_out_discarded = [path_out filesep 'discarded'];
    mkdir(path_out_discarded)
end

close all
fig_dim =1*[20 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
colormap gray
w_plot = 10;


if plot_discarded
    display('Discarding spots...')
    for i=1:N_peaks_raw
        if ~accepted(i) % discarded spot
            message = {'Spotsize red: OK', 'Spotsize green: OK'};
            if criteria(i,1)==0
                message{1} = 'Spotsize red: BAD';
            end
            if criteria(i,2)==0
                message{2} = 'Spotsize green: BAD';
            end
            
            
            x_1 = ch1_fit_raw(i,1);
            y_1 = ch1_fit_raw(i,2);

            x_2 = ch2_fit_raw(i,1);
            y_2 = ch2_fit_raw(i,2);

            subplot(1, 2, 1) %red
            plot_subframe(avg_img{peaks_raw(i, 5), 1}, x_1, y_1, w_plot), hold on
            plot(x_1, y_1, 'rx')
            plot(x_2, y_2, 'g+')
            viscircles(ch1_fit_raw(i,1:2), ch1_fit_raw(i,3),'EdgeColor', 'r'  );
            title({['Pair ' num2str(i) ' of '  num2str(N_peaks_raw) ' at (' num2str(round(x_1)) ',' num2str(round(y_1)) ') in ' channel{1} ' channel'], message{1}})
            axis square
            hold off

            subplot(1, 2, 2) %green
            plot_subframe(avg_img{peaks_raw(i, 5), 2}, x_2, y_2, w_plot), hold on
            plot(x_1, y_1, 'rx')
            plot(x_2, y_2, 'g+')
            %viscircles(gca, ch2_fit_raw(i,1:2), ch2_fit_raw(i,3), 'g'  );
            title({['Pair ' num2str(i) ' of '  num2str(N_peaks_raw) ' at (' num2str(round(x_2)) ',' num2str(round(y_2)) ') in ' channel{2} ' channel'], message{2}})
            axis square
            hold off

            print(cur_fig, '-dtiff', '-r150',  [path_out_discarded filesep 'Discarded_' num2str(i) '.tif'])
        end 
    end
end

display(['Discarded ' num2str(sum(~accepted)) ' spot.'])

close all

N_peaks = size(peaks,1);
%}


%% error analysis:
cur_fig = figure('Visible','off', 'PaperPositionMode', 'manual','PaperUnits','centimeters',...
    'PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
%red channel
xhist = 0:0.05:8;

        subplot(3,2,1)
        n=hist( ch1_2fit_err(:,1),xhist);
        bar(xhist, n);
        legend({'x err'})
        set(gca, 'XLim', [0 0.5])
        
        subplot(3,2,2)
        n=hist( ch1_2fit_err(:,2),xhist);
        bar(xhist, n);
        legend({'y err'})
        set(gca, 'XLim', [0 0.5])
        
        subplot(3,2,3)
        n=hist( ch1_2fit_err(:,3),xhist);
        bar(xhist, n);
        legend({'sigma err'})
        set(gca, 'XLim', [0 0.5])
        
        subplot(3,2,5)
        hist( ch1_2fit_err(:,4));
        legend({'I-bg err'})
        
        subplot(3,2,6)
        hist( ch1_2fit_err(:,5));
        %bar(xhist, n);
        legend({'bg err'})
        %set(gca, 'XLim', [5 10])
        
       
        print(cur_fig, '-dpng', '-r96',  [path_out filesep 'Errors_fit_red_ch.png'])

        [err_fitr(1,1), err_fitr(1,2)]=normfit(ch1_2fit_err(:,1));
        [err_fitr(2,1), err_fitr(2,2)]=normfit(ch1_2fit_err(:,2));
        [err_fitr(3,1), err_fitr(3,2)]=normfit(ch1_2fit_err(:,3));

%green channel
xhist = 0:0.05:8;

        subplot(3,2,1)
        n=hist( ch2_2fit_err(:,1),xhist);
        bar(xhist, n);
        legend({'x err'})
        set(gca, 'XLim', [0 0.5])
        
        subplot(3,2,2)
        n=hist( ch2_2fit_err(:,2),xhist);
        bar(xhist, n);
        legend({'y err'})
        set(gca, 'XLim', [0 0.5])
        
        subplot(3,2,3)
        n=hist( ch2_2fit_err(:,3),xhist);
        bar(xhist, n);
        legend({'sigma err'})
        set(gca, 'XLim', [0 0.5])
        
        subplot(3,2,5)
        hist( ch2_2fit_err(:,4));
        legend({'I-bg err'})
        
        subplot(3,2,6)
        hist( ch2_2fit_err(:,5));
        %bar(xhist, n);
        legend({'bg err'})
        %set(gca, 'XLim', [5 10])
        
       
        print(cur_fig, '-dpng', '-r96',  [path_out filesep 'Errors_fit_green_ch.png'])

        [err_fitg(1,1), err_fitg(1,2)]=normfit(ch2_2fit_err(:,1));
        [err_fitg(2,1), err_fitg(2,2)]=normfit(ch2_2fit_err(:,2));
        [err_fitg(3,1), err_fitg(3,2)]=normfit(ch2_2fit_err(:,3));        
        
        
err_fitr;
err_fitg;


%% random selection of subsets of spots
%{
subset = round(rand(N_peaks,1)*0.75);
%subset = ones(N_peaks,1);
xy1_rand = peaks(subset==1,1:2);
xy2_rand = peaks(subset==1,3:4);

N_sel=length(xy1_rand)


%% TRE evaluation for different NN

N_par= 50; % N_peaks

xy_ch1=zeros(N_sel,4);
xy_ch2=zeros(N_sel,4);
%av_red=zeros(N_par-6,2);
%av_green=zeros(N_par-6,2);
%fitr=zeros(N_par-6,2);
%fitg=zeros(N_par-6,2);
d_red=zeros(N_sel,1);
d_green=zeros(N_sel,1);
d2_red=zeros(N_sel,1);
d2_green=zeros(N_sel,1);

tre_r_all=size(N_par-6,1);
tre_g_all=size(N_par-6,1);

h = waitbar(0,'Calculating mapping function... please wait');


for i=6:N_par
    for k=1:N_sel 
        xy1_tmp=xy1_rand;
        xy1_tmp(k,:)=[];
        
        xy2_tmp=xy2_rand;
        xy2_tmp(k,:)=[];

            tform_2TO1 = fitgeotrans(xy2_tmp, xy1_tmp,'lwm', i); %moving_points, fixed_points, local weighted mean, nearest neighbours, 1 = f(2)
            tform_1TO2 = fitgeotrans(xy1_tmp, xy2_tmp, 'lwm', i); %moving_points, fixed_points, local weighted mean, nearest neighbours, 2 = f(1)

            % Mapp coordinates
            xy_tmp = transformPointsInverse(tform_1TO2, xy2_rand(k,:));  %%this takes coords in ch2 and transfroms them to coords in ch1
            xy_ch1(k,:) = [xy1_rand(k,:) xy_tmp ];

            xy_tmp = transformPointsInverse(tform_2TO1, xy1_rand(k,:));  %%this takes coords in ch1 and transfroms them to coords in ch2
            xy_ch2(k,:) = [xy_tmp xy2_rand(k,:) ];

            waitbar( (i-5)/(N_par-6) , h, ['Calculating mapping with ' num2str(i) ' of ' num2str(N_par) ' nn for spot ' num2str(k) ' of ' num2str(N_sel) ' done']) % update waitbar

    end


d2_red(:) = (xy_ch1(:,1)-xy_ch1(:,3)).^2 + (xy_ch1(:,2)-xy_ch1(:,4)).^2; % dist on red ch 
d2_green(:) = (xy_ch2(:,1)-xy_ch2(:,3)).^2 + (xy_ch2(:,2)-xy_ch2(:,4)).^2; % dist on green ch
    
d_red(:) = sqrt(d2_red); % dist on red ch
d_green(:) = sqrt(d2_green); % dist on green ch

tre_r_all(i-5)=sqrt((sum(d2_red))/(N_sel));
tre_g_all(i-5)=sqrt((sum(d2_green))/(N_sel));

%av_red(i,1)=mean(d_red);
%av_red(i,2)=median(d_red);

%av_green(i,1)=mean(d_green);
%av_green(i,2)=median(d_green);

%[fitr(i-5,1), fitr(i-5,2)]=normfit(d_red);
%[fitg(i-5,1), fitg(i-5,2)]=normfit(d_green);
 
end 


%% plot TRE
 nn= [6:1:50];
 nn2= [6:0.1:50];
 y=0.1;
 
subplot(3,1,1)
plot(nn, tre_r_all, 'rx', 'MarkerSize',6)
hold on 
plot(nn2,y, 'b- ','linewidth',1)
set(gca, 'ylim', [0.05 0.3])
%set(gca, 'xlim', [8 50])
xlabel('Number of nearest neighbor points');
ylabel('TRE value');

subplot(3,1,2)
plot(nn, tre_g_all, 'go', 'MarkerSize',6)
hold on 
plot(nn2,y, 'b- ','linewidth',1)
set(gca, 'ylim', [0.05 0.3])
%set(gca, 'xlim', [8 50])
xlabel('Number of nearest neighbor points');
ylabel('TRE value');

subplot(3,1,3)
plot(nn, tre_r_all, 'rx', 'MarkerSize',6)
hold on
plot(nn, tre_g_all, 'go', 'MarkerSize',6)
set(gca, 'ylim', [0.05 0.3])
%set(gca, 'xlim', [8 50])
xlabel('Number of nearest neighbor points');
ylabel('TRE value');


print(cur_fig, '-dtiff', '-r150',  [path_out filesep 'NN_test_6to50_' num2str(N_sel) '_zoom.tif'])


%% plot all spots used for TRE
close all
plot(xy1_rand(:,1),xy1_rand(:,2),'rx')
set(gca, 'ylim', [0 512])
set(gca, 'xlim', [0 512])
print(cur_fig, '-dtiff', '-r150',  [path_out filesep 'NN_test_6to50_' num2str(N_sel) '_spots.tif'])

%}
%% Mapping using fitgeotrans
peaks_in_range = peaks;
for i = 1:4
    peaks_in_range(peaks_in_range(:,i)<1 | peaks_in_range(:,i)>512,:) = [];
end
N_peaks = size(peaks_in_range,1);
%subset = round(rand(N_peaks,1)*1);
subset = ones(N_peaks,1);

xy_1 = peaks_in_range(subset==1,1:2); % red coord
xy_2 = peaks_in_range(subset==1,3:4);% green coord
n_mov = peaks_in_range(subset==1,5);
N_sel=length(xy_1);

%%
nn=40;

tform_2TO1 = fitgeotrans(xy_2,xy_1,'lwm', nn); %moving_points, fixed_points, local weihted mean, nearest naeighbours, 1 = f(2)
tform_1TO2 = fitgeotrans(xy_1,xy_2,'lwm', nn); %moving_points, fixed_points, local weihted mean, nearest naeighbours, 2 = f(1)

% Mapp coordinates
xy_tmp = transformPointsInverse(tform_2TO1, xy_1);  %this takes coords in ch1 and transfroms them to coords in ch2
xy_ch2 = [xy_tmp xy_2 n_mov];

%sum(sum((xy_1-xy_2).^2))
%sum(sum(( transformPointsInverse(tform_2TO1, xy_1) -xy_2).^2))


xy_tmp = transformPointsInverse(tform_1TO2, xy_2);  %%this takes coords in ch2 and transfroms them to coords in ch1
xy_ch1 = [xy_1 xy_tmp n_mov];

xy_raw = [xy_1 xy_2 n_mov];
%N_pairs = size(xy_ch1,1);


%% SAVE everything
save([path_out filesep 'data.mat'])
tform = tform_2TO1;
save([path_out filesep 'tform_ ' channel{2} '2' channel{1} '_' num2str(nn) '_nn used_' num2str(N_sel) '_points used .mat'], 'tform')
tform = tform_1TO2;
save([path_out filesep 'tform_ ' channel{1} '2' channel{2} '_' num2str(nn) '_nn used_' num2str(N_sel) '_points used .mat'], 'tform')

%% Plot all spots
close all
cur_fig = figure('Visible','off', 'PaperPositionMode', 'manual','PaperUnits','centimeters',...
    'PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
    

    subplot(1, 2, 1)  %red ch
    plot(xy_ch1(:,1), xy_ch1(:,2),'r+', 'MarkerSize', 5 )
    set(gca, 'XLim', [-20 size(avg_img{1,1},1)+20])
    set(gca, 'YLim', [-20 size(avg_img{1,1},1)+20])
    axis image
    title([num2str(N_sel) 'Beads used for mapping, red channel'])
    
    
    subplot(1, 2, 2)  %green ch
    plot(xy_ch2(:,1), xy_ch2(:,2),'g+', 'MarkerSize', 5 )
    set(gca, 'XLim', [-20 size(avg_img{1,1},1)+20])
    set(gca, 'YLim', [-20 size(avg_img{1,1},1)+20])
    axis image
    title([num2str(N_sel) 'Beads used for mapping, green channel'])
    
    print(cur_fig, '-dpng', '-r96', [path_out filesep 'N spots' num2str(N_sel) ' used for mapping.png'])
%% TRE evaluation for Mapping 

xy_ch1_tre=zeros(N_sel,4);
xy_ch2_tre=zeros(N_sel,4);
d_red=zeros(N_sel,1);
d_green=zeros(N_sel,1);
d2_red=zeros(N_sel,1);
d2_green=zeros(N_sel,1);

h = waitbar(0,'Calculating TRE values... please wait');

for k=1:N_sel 
    xy1_tmp=xy_1;
    xy1_tmp(k,:)=[];

    xy2_tmp=xy_2;
    xy2_tmp(k,:)=[];

    tform_2TO1 = fitgeotrans(xy2_tmp, xy1_tmp,'lwm', nn); %moving_points, fixed_points, local weighted mean, nearest neighbours, 1 = f(2)
    tform_1TO2 = fitgeotrans(xy1_tmp, xy2_tmp, 'lwm', nn); %moving_points, fixed_points, local weighted mean, nearest neighbours, 2 = f(1)

    % Mapp coordinates
    xy_tmp = transformPointsInverse(tform_1TO2, xy_2(k,:));  %%this takes coords in ch2 and transfroms them to coords in ch1
    xy_ch1_tre(k,:) = [xy_1(k,:) xy_tmp ];

    xy_tmp = transformPointsInverse(tform_2TO1, xy_1(k,:));  %%this takes coords in ch1 and transfroms them to coords in ch2
    xy_ch2_tre(k,:) = [xy_tmp xy_2(k,:) ];
    
     waitbar( (k)/(N_sel) , h, ['Calculating mapping with ' num2str(nn) ' nn for spot ' num2str(k) ' of ' num2str(N_sel) ' done']) % update waitbar

    
end
    



d2_red(:) = (xy_ch1_tre(:,1)-xy_ch1_tre(:,3)).^2 + (xy_ch1_tre(:,2)-xy_ch1_tre(:,4)).^2; % dist on red ch 
d2_green(:) = (xy_ch2_tre(:,1)-xy_ch2_tre(:,3)).^2 + (xy_ch2_tre(:,2)-xy_ch2_tre(:,4)).^2; % dist on green ch
    
d_red(:) = sqrt(d2_red); % dist on red ch
d_green(:) = sqrt(d2_green); % dist on green ch

tre_r=sqrt((sum(d2_red))/(N_sel));
tre_g=sqrt((sum(d2_green))/(N_sel));


% write info on a file #2
res = [N_sel, tre_r, tre_g, nn];

file = fopen([ path_out filesep ' res_bs.txt'], 'a' );
  fprintf(file,'N_peaks, tre_r, tre_g,nn, \n');
  fprintf( file, '%6.0f %6.3f %6.3f %6.0f', res);
fclose(file);



%% plot pairs
close all
w_plot = 15;

path_out_pairs = [path_out filesep 'Beads pairs_' num2str(w_fit) 'pxl_w_fit'];
mkdir(path_out_pairs)

h = waitbar(0,'Printing spots... please wait');

fig_dim =[20 10];
cur_fig = figure('Visible','off', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

for i=1:N_sel
    
    
        x1=xy_ch1(i,1);
        y1=xy_ch1(i,2);
        
        x2m=xy_ch1(i,3);
        y2m=xy_ch1(i,4);
        
        x1m=xy_ch2(i,1);
        y1m=xy_ch2(i,2);
        
        x2=xy_ch2(i,3);
        y2=xy_ch2(i,4);
    

            subplot(1, 2, 1) %image in channel 1
            
            plot_subframe(avg_img{peaks(i,5),1}, x1, y1, w_plot), hold on
            plot(x1, y1, 'rx',  x2m,y2m, 'gx',  'MarkerSize', 15 )
            rectangle('Position', [ x1-w_fit-0.5,  y1-w_fit-0.5, 2*w_fit+1, 2*w_fit+1], 'EdgeColor', 'r')
            
            title(['Pair ' num2str(i) ' of '  num2str(N_pairs) ' at (' num2str(round(x1)) ',' num2str(round(y1)) ') in ' channel{1} ' ch'])
            % legend({channel{1} , [channel{2} ' mapped'], [ channel{2} ' unmapped']})
            axis square
            hold off


            subplot(1, 2, 2) %image in channel 2

            plot_subframe(avg_img{peaks(i,5), 2}, x2, y2, w_plot), hold on
            plot(x1m, y1m, 'rx', x2, y2, 'gx', 'MarkerSize', 15 )
            rectangle('Position', [ x2-w_fit-0.5,  y2-w_fit-0.5, 2*w_fit+1, 2*w_fit+1], 'EdgeColor', 'g')

            title(['Pair ' num2str(i) ' of '  num2str(N_pairs) ' at (' num2str(round(x2)) ',' num2str(round(y2)) ') in ' channel{2} ' ch'])
            % legend({channel{1} , [channel{2} ' mapped'], [ channel{2} ' unmapped']})
            axis square
            hold off


            print(cur_fig, '-dtiff', '-r150',  [path_out_pairs filesep 'Pair_' num2str(i) '.tif'])

            waitbar( i/N_pairs , h, ['Printing spot... ' num2str(i) ' of ' num2str(N_pairs) ' done']) % update waitbar

         

end             
display('done printing pairs')

close all


%% Write images
path_out_img = [path_out filesep 'rgb_images'];
mkdir(path_out_img)

for i=1:N_movie
    rgb = zeros(size(avg_img{i,1},1), size(avg_img{i,1},1), 3, 'uint16');
        
    A = min_max_uint16(avg_img{i,1});
    B = min_max_uint16(avg_img{i,2});
    B_tformed= imwarp(B, tform_2TO1, 'OutputView', imref2d(size(A)));
    A_tformed= imwarp(A, tform_1TO2, 'OutputView', imref2d(size(B)));

   rgb(:,:,1) = A;
   rgb(:,:,2) = B_tformed;    
   imwrite(rgb, [path_out_img filesep 'image_' sprintf('%02i', i) '_on_' channel{1} '.tif'])
   
   
   rgb(:,:,1) = A_tformed;
   rgb(:,:,2) = B;    
   imwrite(rgb, [path_out_img filesep 'image_' sprintf('%02i', i) '_on_' channel{2} '.tif'] )
end

%%
% summed images
A = zeros(size(avg_img{1,1}));
B = zeros(size(avg_img{1,1}));
for i=1:N_movie     
    A = A + avg_img{i,1};
    B = B + avg_img{i,2};
end
A = min_max_uint16(A);
B = min_max_uint16(B);

B_tformed= imwarp(B, tform_2TO1, 'OutputView', imref2d(size(A)));
A_tformed= imwarp(A, tform_1TO2, 'OutputView', imref2d(size(B)));
rgb = zeros(size(avg_img{i,1},1), size(avg_img{i,1},1), 3, 'uint16');
rgb(:,:,1) = A;
rgb(:,:,2) = B_tformed;    
imwrite(rgb, [path_out_img filesep 'image_sum_' sprintf('%02i', i) '_on_' channel{1} '.tif'])


rgb(:,:,1) = A_tformed;
rgb(:,:,2) = B;    
imwrite(rgb, [path_out_img filesep 'image_sum_' sprintf('%02i', i) '_on_' channel{2} '.tif'] )

rgb(:,:,1) = A;
rgb(:,:,2) = B;    
imwrite(rgb, [path_out_img filesep 'image_sum_unmapped.tif'] )



%% PLOT distances

display('------------------------- PLOTTING DATA --------------------------')

% plot distances before and after mapping

xhist = 0:0.1:8;
close all
d=zeros(N_pairs,3);

index = 1:N_pairs;
d(:,1) = sqrt((xy_ch1(:,1)-xy_ch2(:,3)).^2 + (xy_ch1(:,2)-xy_ch2(:,4)).^2); % dist unmapped
d(:,2) = sqrt((xy_ch1(:,1)-xy_ch1(:,3)).^2 + (xy_ch1(:,2)-xy_ch1(:,4)).^2); % dist on red ch
d(:,3) = sqrt((xy_ch2(:,1)-xy_ch2(:,3)).^2 + (xy_ch2(:,2)-xy_ch2(:,4)).^2); % dist on green ch

subplot(3, 1, 1)
n = hist(d(:,1), xhist);
bar(xhist, n)
title(['Distances in pxl UNMAPPED'])
set(gca, 'xlim', [0 2])

subplot(3, 1, 2)
n = hist(d(:,2), xhist);
bar(xhist, n)
title(['Distances in pxl, mapping on red ch'])
set(gca, 'xlim', [0 2])

subplot(3, 1, 3)
n = hist(d(:,3), xhist);
bar(xhist, n)
title(['Distances in pxl, mapping on green ch'])
set(gca, 'xlim', [0 2])

print(cur_fig, '-dtiff', '-r300', [path_out filesep 'Mapped_distances.tif'])


%% Plot all spots

    close all

    subplot(1, 2, 1)  %red ch
    plot(xy_ch1(:,1), xy_ch1(:,2),'r+', 'MarkerSize', 5 )
    set(gca, 'XLim', [-20 size(avg_img{1,1},1)+20])
    set(gca, 'YLim', [-20 size(avg_img{1,1},1)+20])
    axis image
    title(['Fiducial used for mapping, red channel'])
    
    
    subplot(1, 2, 2)  %green ch
    plot(xy_ch2(:,1), xy_ch2(:,2),'g+', 'MarkerSize', 5 )
    set(gca, 'XLim', [-20 size(avg_img{1,1},1)+20])
    set(gca, 'YLim', [-20 size(avg_img{1,1},1)+20])
    axis image
    title(['Fiducial used for mapping, green channel'])
    
    print(cur_fig, '-dtiff', '-r300', [path_out filesep 'All spots.tif'])
%%

%{
%% plot each spot
close all
button = questdlg('Plot each pair?','Plot pair','Yes','No','No');
w_plot = 10;
if strcmp(button, 'Yes')
    path_out_pairs = [path_out filesep 'pairs'];
    mkdir(path_out_pairs)
    
    h = waitbar(0,'Printing spots... please wait');

    fig_dim =[20 10];
    cur_fig = figure('Visible','off', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
    for i=1:N_peaks
 
                subplot(1, 2, 1) %image in channel 1
                x = xy_ch1(i,1);
                y = xy_ch1(i,2);
                plot_subframe(images{xy_ch1(i,5), 1}, x, y, w_plot), hold on
                plot(x, y, 'rx', xy_ch1(i,3), xy_ch1(i,4), 'g+', xy_raw(i,3), xy_raw(i,4), 'g.', 'MarkerSize', 15 )
                title(['Pair ' num2str(i) ' of '  num2str(N_peaks) ' at (' num2str(round(x)) ',' num2str(round(y)) ') in ' channel{1} ' channel'])
                legend({channel{1} , [channel{2} ' mapped'], [ channel{2} ' unmapped']})
                axis square
                hold off

                subplot(1, 2, 2) %image in channel 2
                x = xy_ch2(i,3);
                y = xy_ch2(i,4);
                plot_subframe(images{xy_ch2(i,5), 2}, x, y, w_plot), hold on
                plot(x, y, 'gx', xy_ch2(i,1), xy_ch2(i,2), 'r+', xy_raw(i,1), xy_raw(i,2), 'r.', 'MarkerSize', 15 )
                title(['Pair ' num2str(i) ' of '  num2str(N_peaks) ' at (' num2str(round(x)) ',' num2str(round(y)) ') in ' channel{2} ' channel'])
                legend({channel{2} , [channel{1} ' mapped'], [ channel{1} ' unmapped']})
                axis square
                hold off

               
                print(cur_fig, '-dtiff', '-r150',  [path_out_pairs filesep 'Pair_' num2str(i) '.tif'])
                
                waitbar( i/N_peaks , h, ['Printing spot... ' num2str(i) ' of ' num2str(N_peaks) ' done']) % update waitbar


    end             
    display('done printing pairs')
end
close(h)



%% print distance scatterplot and get width
close all
fig_dim =[20 20];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

dx = xy_ch1(:,3)-xy_ch1(:, 1);
dy = xy_ch1(:,4)-xy_ch1(:, 2);
lim = 1.1*[min([min(dx) min(dy)]) max([max(dx) max(dy)]) ];

xhist = lim(1):(lim(2)-lim(1))/10:lim(2);

n_dx = hist(dx, xhist);
n_dy = hist(dy, xhist);


subplot(2, 2, 1)
plot(dx, dy, 'r.', 'MarkerSize', 5)
set(gca, 'XLim', lim)
set(gca, 'YLim', lim)
axis square
ylabel('Distance [pixel]')

subplot(2, 2, 3)
bar(xhist, n_dx, 'r'), hold on
set(gca, 'XLim', lim)
set(gca, 'Ydir', 'reverse')
axis square
legend(['sigma_x= ' num2str(std(dx)) ])
xlabel('Distance [pixel]')

subplot(2, 2, 2)
barh(xhist, n_dy, 'r'), hold on
legend(['sigma_y= ' num2str(std(dy)) ])
set(gca, 'YLim', lim)
axis square

print(cur_fig, '-dtiff', '-r300', [path_out filesep 'Distance_scatter_' channel{1} '.tif'])






%% distance distribution
close all

d = sqrt((xy_raw(:,1)-xy_raw(:,3)).^2 + (xy_raw(:,2)-xy_raw(:,4)).^2);
d_tform = sqrt(   (xy_ch1(:,1)-xy_ch1(:,3)).^2 + (xy_ch1(:,2)-xy_ch1(:,4)).^2   );


xhist = 0:max(d)/20:max(d);
n = hist(d, xhist);

xhist2 = 0:max(d)/200:max(d);
n_tform = hist(d_tform, xhist2);

fig_dim =[20 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

subplot(2, 1, 1)
bar(xhist, n)
set(gca, 'XLim', [xhist(1) xhist(end)])
ylabel('Frequency')
legend({'no correction'})
title('Channel 1')


subplot(2, 1, 2)
bar(xhist2, n_tform)
set(gca, 'XLim', [xhist(1) xhist(end)])
xlabel('Distance [pixel]')
ylabel('Frequency')
legend({'mapped'})
print(cur_fig, '-dtiff', '-r300', [path_out filesep 'Distance_distribution_' channel{1} '.tif'])


%% plot histogram of width sigma 

fig_dim =[15 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

xhist = 0:0.05:4;
[n_ch1] = hist(sqrt(ch1_fit(:,3).^2 + ch1_fit(:,4).^2),xhist);
[n_ch2] = hist(sqrt(ch2_fit(:,3).^2 + ch2_fit(:,4).^2),xhist);
subplot(3,1, 1)
plot(xhist, n_ch1, '-r', xhist, n_ch2, '-g', 'MarkerSize', 15)
set(gca, 'XLim', [xhist(1) xhist(end)])
title('Width of peak')
ylabel('Frequency')
legend({[channel{1} ' channel'], [ channel{2} ' channel']})

subplot(3,1, 2)
bar(xhist, n_ch1, 'r')
set(gca, 'XLim', [xhist(1) xhist(end)])
ylabel('Frequency')
legend({[channel{1} ' channel']})

subplot(3, 1, 3)
bar(xhist, n_ch2, 'g')
set(gca, 'XLim', [xhist(1) xhist(end)])
ylabel('Frequency')
legend({[channel{2} ' channel']})

ylabel('Frequency')
xlabel('$\sqrt{\sigma_x^2 + \sigma_y^2}$ in pixel', 'Interpreter', 'LaTex')

print(cur_fig, '-dtiff', '-r300',  [path_out filesep 'width.tif'])



%% scatterplot of width
close all
fig_dim =[10 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
plot(ch1_fit(:,3), ch1_fit(:,4), 'r.', ch2_fit(:,3), ch2_fit(:,4), 'g.', 'MarkerSize', 15)
xlabel('$\sigma_x$ in pixel', 'Interpreter', 'LaTex')
ylabel('$\sigma_y$ in pixel', 'Interpreter', 'LaTex')
title('Width - scatterplot')
legend({[channel{1} ' channel'], [channel{2} ' channel']})

lim = [0.9*min(min([ch1_fit(:,3) ch1_fit(:,4) ch2_fit(:,3) ch2_fit(:,4)])) 1.1*max(max([ch1_fit(:,3) ch1_fit(:,4) ch2_fit(:,3) ch2_fit(:,4)]))];


set(gca, 'Xlim', lim)
set(gca, 'Ylim', lim)
axis square
print(cur_fig, '-dtiff', '-r300',   [path_out filesep 'width_scatterplot.tif'])




%% print a vector-plot , unscaled
fig_dim =[30 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

subplot(1, 2, 1)
dx = xy_ch2(:,1)-xy_ch1(:, 1);
dy = xy_ch2(:,2)-xy_ch1(:, 2);
quiver(xy_ch1(:,1), xy_ch1(:,2), dx, dy, channel{1}(1) ,'AutoScale','off')
set(gca, 'XLim', [-20 size(images{1,1},1)+20])
set(gca, 'YLim', [-20 size(images{1,1},1)+20])
axis image
title(['Transformation of ' channel{1} ' positions'])

subplot(1, 2, 2)
dx = xy_ch1(:,3)-xy_ch2(:, 3);
dy = xy_ch1(:,4)-xy_ch2(:, 4);
quiver(xy_ch2(:,3), xy_ch2(:,4), dx, dy, channel{2}(1) ,'AutoScale','off')
set(gca, 'XLim', [-20 size(images{1,1},1)+20])
set(gca, 'YLim', [-20 size(images{1,1},1)+20])
axis image
title(['Transformation of ' channel{2} ' positions'])

print(cur_fig, '-depsc2', [path_out filesep 'Transformation_Image_absolute.eps'])

%% print a vector-plot , scaled
fig_dim =[30 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

subplot(1, 2, 1)
dx = xy_ch2(:,1)-xy_ch1(:, 1);
dy = xy_ch2(:,2)-xy_ch1(:, 2);
quiver(xy_ch1(:,1), xy_ch1(:,2), dx, dy, channel{1}(1) ,'AutoScale','on')
set(gca, 'XLim', [-20 size(images{1,1},1)+20])
set(gca, 'YLim', [-20 size(images{1,1},1)+20])
axis image
title(['Transformation of ' channel{1} ' positions'])

subplot(1, 2, 2)
dx = xy_ch1(:,3)-xy_ch2(:, 3);
dy = xy_ch1(:,4)-xy_ch2(:, 4);
quiver(xy_ch2(:,3), xy_ch2(:,4), dx, dy, channel{2}(1) ,'AutoScale','on')
set(gca, 'XLim', [-20 size(images{1,1},1)+20])
set(gca, 'YLim', [-20 size(images{1,1},1)+20])
axis image
title(['Transformation of ' channel{2} ' positions'])

print(cur_fig, '-depsc2', [path_out filesep 'Transformation_Image_scaled.eps'])


%%
[X,Y] = meshgrid(50:10:450,50:10:450);
xy = [X(:) Y(:)];


close all
fig_dim =[30 30];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

xy_map = transformPointsInverse(tform_2TO1, xy);  %%this takes coords in ch1 and transfroms them to coords in ch2

dx = xy_map(:,1)-xy(:, 1);
dy = xy_map(:,2)-xy(:, 2);
quiver(xy(:,1), xy(:,2), dx, dy, channel{1}(1) ,'AutoScale','off', 'LineWidth', 1)
set(gca, 'XLim', [-20 size(images{1,1},1)+20])
set(gca, 'YLim', [-20 size(images{1,1},1)+20])
axis image
title(['Transformation of ' channel{1} ' to ' channel{2}])

print(cur_fig, '-depsc2', [path_out filesep 'Transformation_Image_grid.eps'])
%}

%close all
display('Beads mappping finished.')

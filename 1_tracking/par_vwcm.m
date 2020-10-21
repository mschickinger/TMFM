function [ data, vwcm_output ] = par_vwcm(data_path)
% Parallel particle tracking function for two-channel movie objects
% VWCM position estimation + first processing steps (mapping,
% averaging, rms filtering...)
% result is stored in 'data' struct and saved to disk

%% Load data
tic
cd(data_path)
load('data_proc.mat', 'pos_in_frame', 'mapping');
mapping = mapping;
if mapping == 1
    load('data_proc.mat', 'tform');
else
    tform = cell(2,1);
end
load('data_spot_pairs.mat', 'data');
load('movie_objects.mat', 'ch1', 'ch2');
movies = [ch1 ch2];
N_movie = size(movies,1);

%% Create log file
FID = fopen([data_path filesep 'par_vwcm_log.txt'], 'w');
log_string = 'Data loaded. Starting VWCM position estimation. Time elapsed is ';
fprintf(FID, [datestr(now, 'yyyy-mm-dd, HH:MM') ', status:' '\n']);
fprintf(FID, [log_string datestr(toc/86400, 'HH:MM:SS.') '\n']);
fclose(FID);

%% VWCM estimator

log_string = 'VWCM in channel %1d, movie %1d of %1d done. Time elapsed is ';

% Prepare VWCM output array
vwcm_output = cell(size(movies,1),1);
for m=1:size(movies,1)
    vwcm_output{m} = cell(size(pos_in_frame{m,1}{1},1),2);
    for ch = 1:2    
        for s = 1:size(pos_in_frame{m,ch}{1},1)
            vwcm_output{m}{s,ch} = zeros(length(movies{m,ch}.frames),8);
        end
    end
end

% Perform parallel VWCM calculation and reorder
for m=1:size(movies,1)
    for ch = 1:2
        tmp_obj = movies{m,ch};
        tmp_pos = pos_in_frame{m,ch};
        tmp = cell(length(tmp_obj.frames),1);
        parfor p = 1:length(tmp_obj.frames)
            [pos, delta, N, pos_max, v_max, stdevs] = ...
                tmp_obj.vwcm_in_frame(tmp_obj.frames(p), tmp_pos{p}, 5, 0.01, 1000);
            tmp{p} = [pos delta N pos_max v_max stdevs];
        end
        tmp = [tmp{:}];
        for s = 1:size(tmp,1)
            for i = 1:8
                vwcm_output{m}{s,ch}(:,i) = tmp(s,i:8:end);
            end
        end
        FID = fopen('par_vwcm_log.txt', 'a');
        fprintf(FID, [datestr(now, 'yyyy-mm-dd, HH:MM') ', status:' '\n']);
        fprintf(FID, [log_string datestr(toc/86400, 'HH:MM:SS.') '\n'], ch, m, N_movie);
        fclose(FID);
    end
end
cd(data_path)
save -v7.3 'vwcm_output.mat' 'vwcm_output'

%% Output structure part
vwcm_secondary = cell(size(data));
for m = 1:size(data,1)
    tmp = data{m};
    tmp_s = cell(size(data{m}));
    tmp_vwcm_output = vwcm_output{m};
    parfor s = 1:size(data{m},1)
        for ch = 1:2
            % get data from vwcm estimate
            tmp{s,ch}.vwcm.pos = tmp_vwcm_output{s,ch}(:,1:2);
            tmp_s{s,ch}.vwcm.delta = tmp_vwcm_output{s,ch}(:,3);
            tmp_s{s,ch}.vwcm.N = tmp_vwcm_output{s,ch}(:,4);
            tmp_s{s,ch}.vwcm.pos_max = tmp_vwcm_output{s,ch}(:,5:6);
            tmp_s{s,ch}.vwcm.v_max = tmp_vwcm_output{s,ch}(:,7);
            tmp_s{s,ch}.vwcm.stDev = tmp_vwcm_output{s,ch}(:,8);
            
            % complete datasets
            tmp{s,ch}.vwcm.means100 = running_avg_2d_nnz(tmp{s,ch}.vwcm.pos,100);
            tmp{s,ch}.vwcm.medians101 = medfilt1_trunc_2d(tmp{s,ch}.vwcm.pos,101);

            tmp{s,ch}.vwcm.disp100 = tmp{s,ch}.vwcm.pos - tmp{s,ch}.vwcm.means100;
            tmp{s,ch}.vwcm.dispmed101 = tmp{s,ch}.vwcm.pos - tmp{s,ch}.vwcm.medians101;
            
            tmp{s,ch}.vwcm.r = sqrt(tmp{s,ch}.vwcm.disp100(:,1).^2+tmp{s,ch}.vwcm.disp100(:,2).^2);
            
            tmp{s,ch}.vwcm.rms10 = RMSfilt2d(tmp{s,ch}.vwcm.pos,10); 
            
            % mapping
            if mapping
                tmp{s,ch}.vwcm.pos_map = zeros(size(tmp{s,ch}.vwcm.pos));  
                for i = 1:size(tmp{s,ch}.vwcm.pos_map,1)
                    if sum(tmp{s,ch}.vwcm.pos(i,:)) > 0
                        tmp{s,ch}.vwcm.pos_map(i,:) = transformPointsInverse(tform{ch}, tmp{s,ch}.vwcm.pos(i,:));  %%this takes coords in ch2 and transforms them to coords in ch1
                    end
                end
            end
        end
    end
    data{m} = tmp;
    vwcm_secondary{m} = tmp_s;
end

log_string = 'Data structuring done. Saving data. Time elapsed is ';
FID = fopen('par_vwcm_log.txt', 'a');
fprintf(FID, [datestr(now, 'yyyy-mm-dd, HH:MM') ', status:' '\n']);
fprintf(FID, [log_string datestr(toc/86400, 'HH:MM:SS.')]);
fclose(FID);

%% Save data
cd(data_path)
save -v7.3 'data_spot_pairs.mat' 'data'
save -v7.3 -append 'data_archive.mat' 'vwcm_secondary'
display('Data saved')
end


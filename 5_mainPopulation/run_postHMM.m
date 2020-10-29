%% Prepare input for postHMM: 
inputPostHMM.indices = indicesHMM;
inputPostHMM.XY = xyHMMcorr;
inputPostHMM.state_trajectories = state_trajectories;
inputPostHMM.medI = medI;
inputPostHMM.ranges = zeros(size(inputPostHMM.indices));
inputPostHMM.ex_int = cell(size(inputPostHMM.medI));
for i = setdiff(1:size(inputPostHMM.ranges,1),index_discard)
    inputPostHMM.ranges(i,:) = [arxv{i}.segments(1,1) arxv{i}.segments(end,2)];
    if ~isempty(ignore{inputPostHMM.indices(i,1)})
        inputPostHMM.ex_int{i} = ignore{inputPostHMM.indices(i,1)};
    else
        inputPostHMM.ex_int{i} = zeros(0,2);
    end
end
inputPostHMM.indices(index_discard,:) = [];
inputPostHMM.XY(index_discard) = [];
inputPostHMM.state_trajectories(index_discard) = [];
inputPostHMM.medI(index_discard) = [];
inputPostHMM.ranges(index_discard,:) = [];
inputPostHMM.ex_int(index_discard,:) = [];

%% Post-HMM evaluation + Save post-HMM-data 
% !!!Navigate to appropriate folder before!!!

[outputPostHMM] = postHMM(inputPostHMM);
save dataPostHMM.mat outputPostHMM inputPostHMM

% Export Scatter Stats to Igor
if ~exist('SID','var')
    tmp = inputdlg({'Enter sample ID'});
    SID = tmp{1};
end
StatsForIgor = outputPostHMM.scatterStats;
tmp_remove = find(StatsForIgor(:,5) == 0 | StatsForIgor(:,6) == 0);
StatsForIgor(tmp_remove,:) = [];
for i = 3:4
    StatsForIgor(:,i) = StatsForIgor(:,i)./sqrt(StatsForIgor(:,i+2));
end
stats_to_igor(StatsForIgor(:,1:4), SID)
display('Saved .txt file for Igor scatter plot')
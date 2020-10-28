%{
THINGS YOU NEED TO DO BEFOREHAND:
- set the following parameters for each dataset:
    - jumps: jumpMovs, jumpFrames
    - Dmax (maximum physically stretched out excursion in XY; ctr len)
    - sigManual (optional)
%}

%% Extract whole dataset
indicesG = zeros(size(data{1},1),2);
indicesG(:,1) = 1;
indicesG(:,2) = (1:size(indicesG,1))';

%% No jumps in simulated data
jumpMovs = zeros(1,length(data)); % RESET FOR EVERY DATASET
jumpFrames = cell(size(jumpMovs));
% set intervals to be ignored
ignore = cell(size(jumpMovs));

%% Produce cell array with xy-displacements
xyG = cell(size(indicesG,1),1);

for i = 1:numel(xyG)
    xyG{i} = data{1}{i}.vwcm.pos';
end

%% Find number of traces containing more than a certain percentage of frames above increasing threshold levels
% (exclude data points with unrealistic values from statistic)
Dmax = 4; %INPUT - RESET FOR EVERY DATASET
tol = 0.0001;
threshs = 1:0.05:7;
nPmillAbove = zeros(length(threshs),1);
for i = 1:length(xyG)
    tmp_data = xyG{i};
    tmp_data = tmp_data(:,max(tmp_data,[],1)<=Dmax);
    for j = 1:length(threshs)
        tmpN = sum(double(abs(tmp_data(1,:))>threshs(j) | abs(tmp_data(2,:))>threshs(j)));
        nPmillAbove(j) = nPmillAbove(j) + double(tmpN/length(tmp_data)<=tol);
    end
end

% Determine threshold for 'sensible' displacement values
P = 0.2; %INPUT
THR = threshs(find(nPmillAbove>=P*length(indicesG),1));

% Plot this number against these threshold levels
figure
plot(threshs, nPmillAbove)
hold on
plot([threshs(1) threshs(end)], ceil(numel(xyG)*P*[1 1]), 'r--')
plot([THR THR], [0 length(xyG)],'r--')

%% Pick best-suited interval of each trace for HMM evaluation:
intervalsHMM = zeros(size(indicesG));
xyHMM = cell(size(xyG));
for i = 1:size(intervalsHMM,1)
    tmpM = indicesG(i,1);
    tmpS = indicesG(i,2);
    %tmpI = longest_good_interval(xyG{i},Dmax,THR,'ignore',ignore{tmpM});
    tmpI = [1 length(xyG{i})];
    if ~isempty(tmpI)
        intervalsHMM(i,:) = tmpI;
        xyHMM{i} = xyG{i}(:,tmpI(1):tmpI(2));
    end
end
indicesHMM = indicesG(intervalsHMM(:,1)~=0,:);
xyHMM = xyHMM(intervalsHMM(:,1)~=0);
intervalsHMM(intervalsHMM(:,1)==0,:) = [];

%% Correct XY-trajectories for instrumental vibration with sample average
xyHMMcorr = correctXYensemble(xyHMM,indicesHMM,intervalsHMM);
for i = 1:length(xyHMMcorr)
    if all(xyHMMcorr{i}(1,:)==0) || all(xyHMMcorr{i}(2,:)==0)
        display(['Something is wrong with xyHMMcorr at index ' num2str(i)])
    end
end
% if there's only one spot in a movie, the correction doesn't make sense!

%% Fake medI
% maybe 'simulate' division in intensity intervals later on
medI = cell(xyHMM);
for i = 1:length(medI)
    medI{i} = 9001*ones(size(xyHMM{i},2),1);
end

%% Check distribution of data points in intensity intervals
%iEdges = [6750 7500 8500];
iEdges = [7000 8000 9000];
foo2 = N_below(medI,iEdges);
disp(foo2.N/foo2.N_all)

%% Spot-by-spot HMM analysis (dividing each trace into intensity segments)
models = cell(size(xyHMM));
state_trajectories = cell(size(xyHMM));
arxv = cell(size(xyHMM));
sigManual = [0.4 1.2];
h = waitbar(0,['Spot-by-spot HMM analysis: ' num2str(0) ' of ' num2str(length(xyHMM)) ' done.']);
tic
for i = 1:length(xyHMM)
    [models{i}, state_trajectories{i}, arxv{i}] = mlhmmINTsegments(xyHMMcorr{i}, medI{i}, iEdges, intervalsHMM(i,1), 'sigmas', sigManual);
    waitbar(i/length(xyHMM),h,['Spot-by-spot HMM analysis: ' num2str(i) ' of ' num2str(length(xyHMM)) ' done.']);
end
toc
close(h)

%% Extract suitable starting values sigManual

SIGMA = [];
for i = 1:length(arxv)
    if isfield(arxv{i},'models')
        for j = 1:length(arxv{i}.models)
            for k = 1:2
                SIGMA = [SIGMA arxv{i}.models{j}.states{k}.sigma];
            end
        end
    end
end

figure
histogram(SIGMA,100)

%% Save data from HMM analysis !!!Navigate to appropriate folder before!!!
save HMMdata1.mat state_trajectories arxv iEdges xyHMM xyHMMcorr indicesHMM intervalsHMM medI ignore

%% GUI for inspection of state-assigned trajectories:
%mlmodel = model8_7;
%xyHMM = arxv8_7.xyHMM;
%inDisp = Arxv.indices;
YLIM = [0 3];
%inData = [];
%inData = redo_indices;
%inData = index_trunc;
inData = 1:numel(state_trajectories);
inDisp = indicesHMM(inData,:);
m_start = 1;
i = 1;%find(inDisp(:,1)==m_start,1);
while isempty(state_trajectories{i}) && i<=length(inData)
    i = i+1;
end
ts = figure('Units','normalized','Position',[0 0 1 1]);
for p = 1:4
    subplot(4,1,p)
end
bBack = uicontrol('Style', 'pushbutton', 'String', 'Back', 'Units', 'normalized', 'Position', [0.025 0.8 0.05 0.04], 'Callback', 'if i > 1 i = i-1; end, while isempty(state_trajectories{i}) && i>1 i = i-1; end, uiresume', 'FontSize', 12);
bNext = uicontrol('Style', 'pushbutton', 'String', 'Next', 'Units', 'normalized','Position', [0.925 0.8 0.05 0.04], 'Callback', 'if i < length(inData) i = i+1; end, while isempty(state_trajectories{i}) && i<length(state_trajectories) i = i+1; end, uiresume', 'FontSize', 12);
loLim = uicontrol('Style', 'edit', 'Units', 'normalized', 'Position', [0.025 0.2 0.03 0.03]);
hiLim = uicontrol('Style', 'edit', 'Units', 'normalized', 'Position', [0.06 0.2 0.03 0.03]);
bSet = uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'String', 'Set Xlims', 'Position', [0.025 0.15 0.065 0.04], 'Callback', 'for p = 1:4 subplot(4,1,p), xlim([str2double(loLim.String) str2double(hiLim.String)]); end', 'FontSize', 12);
bReset = uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'String', 'Reset', 'Position', [0.025 0.1 0.065 0.04], 'Callback', 'for p = 1:4 subplot(4,1,p), xlim auto; end', 'FontSize', 12);
bDone = uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'String', 'Done', 'Position', [0.925 0.15 0.05 0.04], 'Callback', 'go_on = 0; uiresume', 'FontSize', 12);

go_on = 1;
while go_on
    tmpXY = arxv{inData(i)}.XY;
    tmpS = state_trajectories{inData(i)};
    tmpRMS = data{inDisp(i,1)}{inDisp(i,2),1}.vwcm.rms10(arxv{inData(i)}.segments(1):arxv{inData(i)}.segments(end));
    plot_twostate(tmpXY,tmpS,tmpRMS');
    subplot(4,1,1)
    hold on
    plot(double(tmpS)+0.75, 'k')
    ylim(YLIM)
    title(['Movie ' num2str(inDisp(i,1)) ', spot ' num2str(inDisp(i,2)), ', index ' num2str(i) '/' num2str(length(inData))],'FontSize',16)
    uiwait(gcf)
end
display('Done.')
close(ts)

%{
%% Consensus model
bad_indices = [25 30 51 53 54 62 83];
consensus_indices = setdiff(1:100,bad_indices);
consensus_models = make_consensus_models(state_trajectories(consensus_indices),xyHMMcorr(consensus_indices), iEdges, medI(consensus_indices));
save consensus_models.mat consensus_models
%% Redo analysis for failed trajectories:
redo_indices = bad_indices;
h = waitbar(0,['Spot-by-spot HMM analysis: ' num2str(0) ' of ' num2str(length(redo_indices)) ' done.']);
tic
for i = 1:length(redo_indices)
    j = redo_indices(i);
    [models{j}, state_trajectories{j}, arxv{j}] = mlhmmINTsegments(xyHMMcorr{j}, medI{j}, iEdges, intervalsHMM(j,1), 'initial_models', consensus_models);
    waitbar(i/length(redo_indices),h,['Spot-by-spot HMM analysis: ' num2str(i) ' of ' num2str(length(redo_indices)) ' done.']);
end
toc
close(h)
%}

%{
%% truncate or discard data from specific particles
% INPUT SPECIFICALLY FOR EVERY NEW DATASET:
discard_manual = [];
truncate_from = ...
    [ ...
];
truncate_to = ...
    [ ...
];
if ~isempty(truncate_from)
    index_truncate_from = truncate_from(:,1); %[,];
    limit_truncate_from = truncate_from(:,2); %[,];
else
    index_truncate_from = [];
    limit_truncate_from = [];
end
if ~isempty(truncate_to)
    index_truncate_to = truncate_to(:,1); %[,];
    limit_truncate_to = truncate_to(:,2); %[,];
else
    index_truncate_to = [];
    limit_truncate_to = [];
end
% adjust HMM intervals, xyHMM, xyHMMcorr, medI
for i = 1:length(index_truncate_from)
    intervalsHMM(index_truncate_from(i),2) = intervalsHMM(index_truncate_from(i),1) + limit_truncate_from(i);
    xyHMM{index_truncate_from(i)} = xyHMM{index_truncate_from(i)}(:,1:limit_truncate_from(i));
    xyHMMcorr{index_truncate_from(i)} = xyHMMcorr{index_truncate_from(i)}(:,1:limit_truncate_from(i));
    medI{index_truncate_from(i)} = medI{index_truncate_from(i)}(1:limit_truncate_from(i));
end
for i = 1:length(index_truncate_to)
    intervalsHMM(index_truncate_to(i),1) = intervalsHMM(index_truncate_to(i),1) + limit_truncate_to(i);
    xyHMM{index_truncate_to(i)} = xyHMM{index_truncate_to(i)}(:,limit_truncate_to(i)+1:end);
    xyHMMcorr{index_truncate_to(i)} = xyHMMcorr{index_truncate_to(i)}(:,limit_truncate_to(i)+1:end);
    medI{index_truncate_to(i)} = medI{index_truncate_to(i)}(limit_truncate_to(i)+1:end);
end


%% Redo analysis in truncated trajectories -> save again
index_trunc = reshape(union(index_truncate_to,index_truncate_from),1,[]);
h = waitbar(0,['Re-doing HMM analysis: ' num2str(0) ' of ' num2str(length(index_trunc)) ' done.']);
tic
for i = index_trunc
    waitbar(find(index_trunc==i)/length(index_trunc),h,['Redoing HMM analysis: Current index is ' num2str(i) ' (' num2str(find(index_trunc==i)) ' of ' num2str(length(index_trunc)) ')']);
    [models{i}, state_trajectories{i}, arxv{i}] = mlhmmINTsegments(xyHMMcorr{i}, medI{i}, iEdges, intervalsHMM(i,1), 'sigmas', sigManual);
end
toc
close(h)
save HMMdata2.mat state_trajectories arxv iEdges xyHMM xyHMMcorr indicesHMM intervalsHMM

discard = zeros(1,length(state_trajectories));
for i = 1:length(discard)
    discard(i) = isempty(state_trajectories{i});
end

index_discard = unique([find(discard==1), discard_manual]);
%save HMMsortout.mat truncate_from truncate_to index_discard
%}
%% Prepare input for postHMM: 
index_discard = [];
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

%% Post-HMM evaluation + save post-HMM-data
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
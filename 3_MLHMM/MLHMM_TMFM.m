%{
THINGS YOU NEED TO DO BEFORE / DURING MLHMM SETUP:
- if the datasets were reduced for long-time storage, run the following
code after loading the data_spot_pairs.mat file
for m = 1:size(data,1)
    for s = 1:size(data{m},1)
        for ch = 1:2
            data{m}{s,ch}.med_itrace = medfilt1_trunc(data{m}{s,ch}.itrace,20);
            data{m}{s,ch}.vwcm.medians101 = medfilt1_trunc_2d(data{m}{s,ch}.vwcm.pos,101);
            data{m}{s,ch}.vwcm.dispmed101 = data{m}{s,ch}.vwcm.pos - data{m}{s,ch}.vwcm.medians101;
            data{m}{s,ch}.vwcm.r = sqrt(data{m}{s,ch}.vwcm.dispmed101(:,1).^2+data{m}{s,ch}.vwcm.dispmed101(:,2).^2);
        end
    end
end
- set the following parameters for each dataset:
    - intervals to be igored: ignore (cell array)
    - Dmax: maximum physically stretched out radial excursion (contour length)
    - iEdges:
        - standard for distal geometries: [7000 8000 9000]
        - standard for proximal geometries: [6500 7500]
    - sigManual (optional): starting values for the x/y standard deviations
%}

%% Set sample ID and load data
clear variables
run('my_prefs.m')

SID = 'C001';
load('data_spot_pairs.mat')
load('GiTSiK.mat')

%% Extract dataset from GiTSiK
tmp = 0;
for m = 1:length(GiTSiK)
   tmp = tmp + sum(GiTSiK.behaviour{m} == 2);
end
indicesG = zeros(tmp,2);
counter = 1;
for m = 1:size(data,1)
    for i = find(GiTSiK.behaviour{m} == 2)';
        indicesG(counter,:) = [m i];
        counter = counter+1;
    end
end

%% Account for jumps in the movies
% Correct the XY-trajectories from movies that have 'jumps' in it
% at every site of a sudden 'jump':
% - rewrite the medians101 arrays around that frame number (+/- 50 frames)
% - recalculate the dispmed101 displacements for the same frame range
% (Truncate traces from first occuring zero onwards)
jumpMovs = zeros(1,length(GiTSiK.behaviour)); % RESET FOR EVERY DATASET
jumpFrames = cell(size(jumpMovs));
% set intervals to be ignored
ignore = cell(size(jumpMovs));
%ignore{2} = [19309 19311; 30050 30100; 30234 30236; 30414 30416]; % RESET FOR EVERY DATASET
if exist('jumps.mat','file')
    display('loading jumps.mat file')
    load jumps.mat
    for m = 1:length(jumpMovs)
        if ~isempty(jumps{m})
            jumpMovs(m) = 1;
            jumpFrames{m} = jumps{m};
            ignore{m} = zeros(numel(jumps{m}),2);
            for j = 1:numel(jumps{m})
                ignore{m}(j,1) = min(jumps{m}(j)-50,floor(jumps{m}(j)/100)*100+1);
                ignore{m}(j,2) = max(jumps{m}(j)+50,ceil(jumps{m}(j)/100)*100);
            end
        end
    end    
end

%% Produce cell array with xy-displacements
xyG = cell(length(indicesG),1);
truncate = []; %[3 82 1 14500]; % RESET FOR EVERY DATASET

for i = 1:length(indicesG)
    if ismember(indicesG(i,1),jumpMovs)
        tmpMF101 = data{indicesG(i,1)}{indicesG(i,2),1}.vwcm.medians101;
        for fNum = reshape(jumpFrames{indicesG(i,1)},1,[])
            for j = 1:2
                tmpMF101(fNum+(-50:-1),j) = tmpMF101(fNum-51,j);
                tmpMF101(fNum+(0:49),j) = tmpMF101(fNum+50,j);
            end
        end
        tmp_data = data{indicesG(i,1)}{indicesG(i,2),1}.vwcm.pos' - tmpMF101';
    else
        tmp_data = data{indicesG(i,1)}{indicesG(i,2),1}.vwcm.dispmed101';
    end
    tmpL = find(data{indicesG(i,1)}{indicesG(i,2),1}.vwcm.r==0,1);
    if isempty(tmpL)
        xyG{i} = tmp_data;
    else
        xyG{i} = tmp_data(:,1:tmpL-1);
    end
end
for j = 1:size(truncate,1)
    tmpI = find(indicesG(:,1)==truncate(j,1) & indicesG(:,2)==truncate(j,2));
    xyG{tmpI} = xyG{tmpI}(:,truncate(j,3):truncate(j,4));
end

%% Find number of traces containing more than a certain percentage of frames above increasing threshold levels
% (exclude data points with unrealistic values from statistic)
Dmax = 7.5; %INPUT - RESET FOR EVERY DATASET
tol = 0.0001;
threshs = 1:0.05:8;
nPmillAbove = zeros(length(threshs),1);
for i = 1:length(xyG)
    tmp_data = xyG{i};
    tmp_data = tmp_data(:,max(tmp_data,[],1)<=Dmax);
    for j = 1:length(threshs)
        tmpN = sum(double(abs(tmp_data(1,:))>threshs(j) | abs(tmp_data(2,:))>threshs(j)));
        nPmillAbove(j) = nPmillAbove(j) + double(tmpN/length(tmp_data)<=tol);
    end
end

%% Determine threshold for 'sensible' displacement values
P = 0.2; %INPUT
THR = threshs(find(nPmillAbove>=P*length(indicesG),1));

% Plot this number against these threshold levels
figure
plot(threshs, nPmillAbove)
hold on
plot([threshs(1) threshs(end)], ceil(length(xyG)*P*[1 1]), 'r--')
plot([THR THR], [0 length(xyG)],'r--')

%% Pick best-suited interval of each trace for HMM evaluation:
intervalsHMM = zeros(size(indicesG));
xyHMM = cell(size(xyG));
for i = 1:length(intervalsHMM)
    tmpM = indicesG(i,1);
    tmpS = indicesG(i,2);
    tmpI = longest_good_interval(xyG{i},Dmax,THR,'ignore',ignore{tmpM});
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

%% cell container for the med_itraces of spots/intervals used for HMM analysis
medI = cell(size(indicesHMM,1),1);
for i = 1:length(medI)
    medI{i} = data{indicesHMM(i,1)}{indicesHMM(i,2),1}.med_itrace(intervalsHMM(i,1):intervalsHMM(i,2));
end

%% Check distribution of data points in intensity intervals
%iEdges = [6750 7500 8500];
iEdges = [7500 8500 10000];
foo2 = N_below(medI,iEdges);
disp(foo2.N/foo2.N_all)

%% Spot-by-spot HMM analysis (dividing each trace into intensity segments)
models = cell(size(xyHMM));
state_trajectories = cell(size(xyHMM));
arxv = cell(size(xyHMM));
sigManual = [0.2 1.15];
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

%% GUI for inspection of state-assigned trajectories:
%mlmodel = model8_7;
%xyHMM = arxv8_7.xyHMM;
%inDisp = Arxv.indices;
YLIM = [0 3];
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
    tmpRMS = RMSfilt2d(tmpXY',11);
%    tmpRMS = data{inDisp(i,1)}{inDisp(i,2),1}.vwcm.rms10(arxv{inData(i)}.segments(1):arxv{inData(i)}.segments(end));
    plot_twostate(tmpXY,tmpS,tmpRMS');
    subplot(4,1,1)
    hold on
    plot(double(tmpS)+.75, 'k')
    ylim(YLIM)
    title(['Movie ' num2str(inDisp(i,1)) ', spot ' num2str(inDisp(i,2)), ', index ' num2str(i) '/' num2str(length(inData))],'FontSize',16)
    uiwait(gcf)
end
display('Done.')
close(ts)

%% Consensus model
%consensus_indices = setdiff(1:50, [16 25 92 11 18 27 31 37 40 49 55 77 107 108 117]);
%consensus_indices = [2,3,5,9,11,19,21,23,25,26,27,28,31,33,35,36,37,39,42,44,45,47,48,49,50];
consensus_indices = 1:4:30;
consensus_models = make_consensus_models(state_trajectories(consensus_indices),xyHMMcorr(consensus_indices), iEdges, medI(consensus_indices));
save consensus_models.mat consensus_models
%% Redo analysis for failed trajectories:
redo_indices = 1:numel(models);
%redo_indices = setdiff(1:numel(state_trajectories),consensus_indices);
%redo_indices = 1:numel(state_trajectories);
h = waitbar(0,['Spot-by-spot HMM analysis: ' num2str(0) ' of ' num2str(length(redo_indices)) ' done.']);
tic
for i = 1:length(redo_indices)
    j = redo_indices(i);
    [models{j}, state_trajectories{j}, arxv{j}] = mlhmmINTsegments(xyHMMcorr{j}, medI{j}, iEdges, intervalsHMM(j,1), 'initial_models', consensus_models);
    waitbar(i/length(redo_indices),h,['Spot-by-spot HMM analysis: ' num2str(i) ' of ' num2str(length(redo_indices)) ' done.']);
end
toc
close(h)

%% Save data from HMM analysis !!!Navigate to appropriate folder before!!!
save HMMdata1.mat state_trajectories arxv iEdges xyHMM xyHMMcorr indicesHMM intervalsHMM medI ignore

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
mode = questdlg('Choose mode', 'Choose mode','regular','consensus','regular'); 
h = waitbar(0,['Re-doing HMM analysis: ' num2str(0) ' of ' num2str(length(index_trunc)) ' done.']);
tic
for i = index_trunc
    waitbar(find(index_trunc==i)/length(index_trunc),h,['Redoing HMM analysis: Current index is ' num2str(i) ' (' num2str(find(index_trunc==i)) ' of ' num2str(length(index_trunc)) ')']);
    switch mode
        case 'regular'
        [models{i}, state_trajectories{i}, arxv{i}] = mlhmmINTsegments(xyHMMcorr{i}, medI{i}, iEdges, intervalsHMM(i,1), 'sigmas', sigManual);
        case 'consensus'
        [models{i}, state_trajectories{i}, arxv{i}] = mlhmmINTsegments(xyHMMcorr{i}, medI{i}, iEdges, intervalsHMM(i,1), 'initial_models', consensus_models);
    end
end
toc
close(h)
save HMMdata2.mat state_trajectories arxv iEdges xyHMM xyHMMcorr indicesHMM intervalsHMM medI

discard = zeros(1,length(state_trajectories));
for i = 1:length(discard)
    discard(i) = isempty(state_trajectories{i});
end
index_discard = unique([find(discard==1), discard_manual]);
save HMMsortout.mat truncate_from truncate_to index_discard

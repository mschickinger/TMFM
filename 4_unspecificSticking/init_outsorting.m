%% LOAD: 
clear variables
run('my_prefs.m')
load('data_spot_pairs.mat')
%load('HMMdata1.mat')
load('HMMdata2.mat', 'arxv', 'iEdges', 'indicesHMM', 'intervalsHMM', ...
        'state_trajectories', 'xyHMMcorr')
load('HMMsortout.mat', 'index_discard')
display('data loaded')  

%% cell container for the med_itraces of spots/intervals used for HMM analysis
medI = cell(size(indicesHMM,1),1);
for i = 1:length(medI)
    medI{i} = data{indicesHMM(i,1)}{indicesHMM(i,2),1}.med_itrace(intervalsHMM(i,1):intervalsHMM(i,2));
end

if ~exist('iEdges','var')
    find_iEdges
else
    display(iEdges,'iEdges')
end

%% Set parameters
%discard_manual = []; % RESET FOR EVERY DATASET
if ~exist('discard_manual','var')
    discard_manual = [];
end
index_discard = union(index_discard,discard_manual);
%iEdges = [8000 7000 6300]; % RESET FOR EVERY DATASET - use result from find_iEdges
RMSintSeg = cell(2,numel(iEdges));
xySeg = cell(2,numel(iEdges));
% distlBelow = cell(2,1);
% distXYBelow = cell(2,1);
% indsBelow = cell(2,1);
% for i = 1:2
%     distXYBelow{i} = zeros(2,0);
%     indsBelow{i} = zeros(0,3);
% end
% statesInds = cell(2,1);
for i = 1:numel(RMSintSeg)
    RMSintSeg{i} = zeros(4,0);
    xySeg{i} = zeros(4,0);
end
%
segments = cell(1,length(medI));
segmInds = cell(1,length(medI));
removInds = cell(2,1);
for k = 1:2
    removInds{k} = cell(length(medI),2);
end
% determine segment edges and indices -> store in cell arrays
for i = 1:length(medI)
    [segments{i}, segmInds{i}] = iSegments(medI{i}, iEdges);
end
% get all sigma values for states in segments -> store in cell array
allSigmas = cell(2,1);
for k = 1:2
    allSigmas{k} = cell(1,numel(iEdges));
    for j = 1:length(allSigmas{k})
        allSigmas{k}{j} = zeros(0,2);
    end
end
for i = setdiff(1:length(segmInds),index_discard)
    for j = 1:length(segmInds{i})
        for k = 1:2
            allSigmas{k}{segmInds{i}(j)} = [allSigmas{k}{segmInds{i}(j)}; arxv{i}.models{j}.states{k}.sigma];
        end
    end
end
% median sigma values for states and segments
medSigmas = cell(size(allSigmas));
for i = 1:2
    medSigmas{i} = zeros(size(allSigmas{i},2),2);
    for j = 1:size(medSigmas{i},1)
        medSigmas{i}(j,:) = median(allSigmas{i}{j});
    end
end
%
% Fill global cell array for segments and states
for i = 1:length(arxv)
    if ~isempty(arxv{i}.segments)
        tmpRMS = data{indicesHMM(i,1)}{indicesHMM(i,2),1}.vwcm.rms10(arxv{i}.segments(1):arxv{i}.segments(end))';
        tmpXY = arxv{i}.XY;
        W = 11;
        tmpXY(1,:) = meanfilt1_trunc(tmpXY(1,:),W);
        tmpXY(2,:) = meanfilt1_trunc(tmpXY(2,:),W);
        if length(tmpRMS)~=length(tmpXY) && length(tmpXY)~=length(state_trajectories{i})
            display([num2str(i) ' of ' num2str(length(medI)) ': ERROR'])
        else
            display([num2str(i) ' of ' num2str(length(medI)) ': OK'])
        end
        if ~isempty(segments{i})
            %tmpS = state_trajectories{i}(1:find(state_trajectories{i}~=state_trajectories{i}(end),1,'last'));
            for k = 1:2
                statesInds{k} = find(state_trajectories{i}==k);
                %statesInds{k} = setdiff(find(state_trajectories{i}==k),removInds{i,2});
            end
            for j = 1:length(segmInds{i})
                for k = 1:2
                    tmpIND = intersect(statesInds{k},segments{i}(j,1):segments{i}(j,2));
                    RMSintSeg{k,segmInds{i}(j)} = [RMSintSeg{k,segmInds{i}(j)} [tmpRMS(tmpIND); ...
                                                                                medI{i}(tmpIND)'; ...
                                                                                i*ones(1,length(tmpIND)); ...
                                                                                tmpIND]];
                    xySeg{k,segmInds{i}(j)} = [xySeg{k,segmInds{i}(j)} [tmpXY(:,tmpIND); ...
                                                                        i*ones(1,length(tmpIND)); ...
                                                                        tmpIND]];  
                end
            end 
        end
    end
end

% Remove implausible data points
RMSmax = 10;
intmin = 6000;
for i = 1:numel(RMSintSeg)
    xySeg{i} = xySeg{i}(:,RMSintSeg{i}(1,:)<RMSmax & RMSintSeg{i}(2,:)>=intmin);
    RMSintSeg{i} = RMSintSeg{i}(:,RMSintSeg{i}(1,:)<RMSmax & RMSintSeg{i}(2,:)>=intmin);  
end
xySegZ = xySeg;
RMSintSegZ = RMSintSeg;

% Get density below 1% threshold for all bound intervals
P = 0.01;
P2 = 0.025;
[globThreshs, gThreshs2] = get_globThreshs(RMSintSeg, P, P2);
densities = cell(length(state_trajectories),2);
noneAboveP2 = cell(size(densities));
%maxlBelow = cell(length(state_trajectories),2);
stateFrames = cell(length(state_trajectories),2);
h = waitbar(0,'');
for isp = setdiff(1:size(densities,1),index_discard)
    waitbar(isp/length(state_trajectories),h,['getting density statistics ' ...
        num2str(isp) ' of ' num2str(length(state_trajectories)) '.']);
    if ~isempty(segments{isp})
        straj = state_trajectories{isp};       
        tmpRMS = data{indicesHMM(isp,1)}{indicesHMM(isp,2),1}.vwcm.rms10(arxv{isp}.segments(1):arxv{isp}.segments(end))';
        tmpXY = arxv{isp}.XY;
        % get starts and lengths of states
        steps = find(diff(straj)~=0) + 1;
        if ~isempty(steps)
            % Start frames and lengths of states
            S = reshape(steps(1:end-1),length(steps)-1,1);
            L = reshape(steps(2:end)-steps(1:end-1),length(steps)-1,1);
            % Assign type of states
            updn = zeros(size(S));
            for i = 1:length(updn)
                updn(i) = sign(straj(steps(i))-straj(steps(i)-1));
            end
            % divide in hi and lo states
            states = cell(1,2);
            states{2} = [S(updn==1) L(updn==1)];
            states{1} = [S(updn==-1) L(updn==-1)];
            for k = 1:2
                states{k}(sum(states{k},2)>segments{isp}(end),:) = [];
            end
            stateFrames(isp,:) = states;
        end
        for k = 1:2
            densities{isp,k} = zeros(size(states{k}));
            noneAboveP2{isp,k} = zeros(size(states{k},1),1);
            tmpSeg = 1;
            for i = 1:size(densities{isp,k},1)
                if ~isnan(states{k}(i,1))
                    while states{k}(i,1) > segments{isp}(tmpSeg,2)
                        tmpSeg = tmpSeg + 1;
                    end
                    tmpA = states{k}(i,2)>1;
                    tmpB = zeros(1,states{k}(i,2));
                    tmpI = states{k}(i,1):sum(states{k}(i,:))-1;
                    tmpF = 0;
                    while sum(states{k}(i,:)) > segments{isp}(tmpSeg,2)
                        tmpA = tmpA && all(tmpRMS(tmpI((tmpF+1):(segments{isp}(tmpSeg,2)-states{k}(i,1))))<gThreshs2(k,segmInds{isp}(tmpSeg)));
                        tmpB((tmpF+1):(segments{isp}(tmpSeg,2)-states{k}(i,1))) = ...
                            tmpRMS(tmpI((tmpF+1):(segments{isp}(tmpSeg,2)-states{k}(i,1))))<=globThreshs(k,segmInds{isp}(tmpSeg));
                        tmpF = segments{isp}(tmpSeg,2)-states{k}(i,1);
                        tmpSeg = tmpSeg + 1;
                    end
                    if tmpF < states{k}(i,2)
                        tmpA = tmpA && all(tmpRMS(tmpI(tmpF+1:end))<gThreshs2(k,segmInds{isp}(tmpSeg)));
                        tmpB(tmpF+1:end) = tmpRMS(tmpI(tmpF+1:end))<=globThreshs(k,segmInds{isp}(tmpSeg));
                    end
                    densities{isp,k}(i,1) = sum(tmpB)/length(tmpB);
                    noneAboveP2{isp,k}(i) = tmpA;
                    Nmax = 100;
                    Lmin = [Nmax, (Nmax+2*floor(W/2))];
                    if states{k}(i,2)>Lmin(k)
                        if k==2
                            tmpB = tmpB(floor(W/2)+1:end-floor(W/2));
                        end
                        tmpB2 = zeros(1,length(tmpB)-Nmax+1);
                        for j = 1:length(tmpB2)
                            tmpB2(j) = sum(tmpB(j:j+Nmax-1));
                        end
                        densities{isp,k}(i,2) = max(tmpB2)/Nmax;
                    elseif k==1
                        densities{isp,k}(i,2) = densities{isp,k}(i,1);
                    else
                        densities{isp,k}(i,2) = NaN;
                    end
                else
                    densities{isp,k}(i,:) = NaN;
                end
            end
        end
    end
end
close(h)
stateFramesZ = stateFrames;
densitiesZ = densities;
noneAboveP2Z = noneAboveP2;
RMSintSegZ = RMSintSeg;
removIndsZ = removInds;
% Check out distribution of densities
[allD, allDmax] = get_allDmax(densities, index_discard);
%
display('Outsorting initialized.')

%%
%
figure('Units','normalized','Position',[0 0 1 1])
for k = 1:2
    subplot(2,2,k)
    histogram(allD{k})
    %ylim([0 20])
    subplot(2,2,k+2)
    histogram(allDmax{k})
    title(['max' num2str(k)])
    %ylim([0 20])
end
%}
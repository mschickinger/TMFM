%% get all sigma values for states -> store in cell array
for k = 2:-1:1
    allSigmas{k} = zeros(0,2);
end
for i = 1:length(arxv)
    for k = 1:2
        allSigmas{k} = [allSigmas{k}; arxv{i}.models{1}.states{k}.sigma];
    end
end
% median sigma values for states and segments
medSigmas = {[0 0],[0 0]};
for k = 1:2
    medSigmas{k} = median(allSigmas{k});
end

%% write all the XY and RMS11 values for the states in cell arrays
W = 11;
RMSall = cell(2,1);
XYall = cell(2,1);
for i = 1:numel(RMSall)
    RMSall{i} = zeros(1,0);
    XYall{i} = zeros(3,0);
end
for i = 1:length(state_trajectories)
    tmpRMS = simRMS11{i};
    tmpXY = simXY{i};
    W = 11;
    tmpXY(1,:) = meanfilt1_trunc(tmpXY(1,:),W);
    tmpXY(2,:) = meanfilt1_trunc(tmpXY(2,:),W);
    if length(tmpRMS)~=length(tmpXY) && length(tmpXY)~=length(state_trajectories{i})
        display([num2str(i) ' of ' num2str(length(simRMS11)) ': ERROR'])
    else
        display([num2str(i) ' of ' num2str(length(simRMS11)) ': OK'])
    end
    for k = 1:2
        tmpIND = find(state_trajectories{i}==k);
        RMSall{k} = [RMSall{k} tmpRMS(tmpIND)];
        XYall{k} = [XYall{k} [tmpXY(:,tmpIND); ...
            sqrt(tmpXY(1,tmpIND).^2 + tmpXY(2,tmpIND).^2).*sqrt(W)./mean(medSigmas{k})]];  
    end
end
%% Remove implausible data points
RMSmax = 10;
for i = 1:numel(RMSall)
    XYall{i} = XYall{i}(:,RMSall{i}(1,:)<RMSmax);
    RMSall{i} = RMSall{i}(:,RMSall{i}(1,:)<RMSmax);  
end

%% Get global thresholds, densities (areas) and distributions

P = 0.01;
globThreshs = get_globThreshs(RMSall, P);
densities = cell(length(state_trajectories),2);
areas = cell(length(state_trajectories),2);
stateFrames = cell(length(state_trajectories),2);
h = waitbar(0,'');
for isp = 1:length(state_trajectories)
    waitbar(isp/length(state_trajectories),h,['getting density and area statistics ' ...
    num2str(isp) ' of ' num2str(length(state_trajectories)) '.']);
    straj = state_trajectories{isp};       
    tmpRMS = simRMS11{isp};
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
        stateFrames(isp,:) = states;
    end
    for k = 1:2
        densities{isp,k} = zeros(size(states{k}));
        areas{isp,k} = zeros(size(states{k}));
        for i = 1:size(densities{isp,k},1)
            tmpI = states{k}(i,1):sum(states{k}(i,:))-1;
            tmpA = tmpRMS(tmpI)-globThreshs(k);
            tmpB = tmpRMS(tmpI)<=globThreshs(k);
            densities{isp,k}(i,1) = sum(tmpB)/length(tmpB);
            areas{isp,k}(i,1) = abs(sum(tmpB.*tmpA));
            Nmax = 100;
            Lmin = [Nmax, (Nmax+2*floor(W/2))];
            if states{k}(i,2)>Lmin(k)
                if k==2
                    tmpB = tmpB(floor(W/2)+1:end-floor(W/2));
                end
                tmpA2 = zeros(1,length(tmpB)-Nmax+1);
                tmpB2 = zeros(1,length(tmpB)-Nmax+1);
                for j = 1:length(tmpB2)
                    tmpA2(j) = abs(sum(tmpA(j:j+Nmax-1).*tmpB(j:j+Nmax-1)));
                    tmpB2(j) = sum(tmpB(j:j+Nmax-1));
                end
                densities{isp,k}(i,2) = max(tmpB2)/Nmax;
                areas{isp,k}(i,2) = max(tmpA2);
            elseif k==1
                densities{isp,k}(i,2) = densities{isp,k}(i,1);
                areas{isp,k}(i,2) = areas{isp,k}(i,1);
            else
                densities{isp,k}(i,2) = NaN;
                areas{isp,k}(i,2) = NaN;
            end
        end
    end
end
close(h)

%%
[allD, allDmax, allA, allAmax] = get_allADmax(densities, [], areas);

limitsD = [0 0];
limitSEM = [0 0]; 
for k = 1:2
    %tmp = sort(allDmax{k});
    limitsD(k) = max(allDmax{k});%tmp(ceil(0.999*length(tmp)));
    tmp = sort(XYall{k}(3,RMSall{k}<globThreshs(k)));
    limitSEM(k) = tmp(ceil(0.999*length(tmp)));
end

save unspecific.mat RMSall XYall densities allD allDmax limitsD limitSEM
display(limitsD,'limitsD')
display(limitSEM,'limitSEM')
%% Plots
figure('Units','normalized','Position',[0 0 1 1])
for k = 1:2
    subplot(2,2,k)
    histogram(allD{k})
    ylim([0 20])
    subplot(2,2,k+2)
    histogram(allDmax{k})
    title(['max' num2str(k)])
    ylim([0 20])
    hold on
    plot(limitsD(k)*[1 1],[0 20],'r--')
end
%%
figure('Units','normalized','Position',[0 0 1 1])
for k = 1:2
    subplot(2,2,k)
    histogram(XYall{k}(3,RMSall{k}<globThreshs(k)))
    subplot(2,2,k+2)
    histogram(XYall{k}(3,RMSall{k}<globThreshs(k)),'normalization','cdf')
    ylim([.99 1])
    hold on
    XLIM = xlim; YLIM = ylim;
    plot(XLIM,.999*[1 1],'r--')
    plot(4*[1 1],YLIM,'r--')
end
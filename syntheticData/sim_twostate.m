clear variables
close all
clear variables
run('my_prefs.m')

%% Generate dummy data
% MAKE SURE YOU ARE INSIDE THE RIGHT FOLDER!!!!!!
SID = 'T000';
sample_means = {50, 50}; % in seconds
Tmin = [1.3 1.0];

k_in = get_corrected_rates(sample_means, Tmin, Inf);
disp(1./k_in)
%%
simParams.tau = [30 20];%1./k_in; %[];
simParams.tpf = 50;
simParams.L = 45000;
simParams.mu = 0;
simParams.sigma = [0.4 1.2];
simParams.wSize = 11;

N = 100;

simStraj = cell(N,1);
simXY = cell(N,1);
simRMS11 = cell(N,1);
simT = cell(N,2);
h = waitbar(0);
rng default %for reproducibility;
for n = 1:N
    waitbar(n/N,h,['creating random trajectory ' ...
        num2str(n) ' of ' num2str(N) '.']);
    [simStraj{n}, simXY{n}, simRMS11{n}, simT(n,:)] = sim_twostate_XY_RMS(simParams);
    %fprintf('%.2f\t%.2f\t%.2f\t%.2f\n',mean(simT{n,1}), mean(simT{n,2}), mean(vertcat(simT{1:n,1})), mean(vertcat(simT{1:n,2})))
end
close(h)
%
data = {cell(N,1)};
for s = 1:N
    data{1}{s}.vwcm.pos = simXY{s}';
    data{1}{s}.vwcm.rms10 = simRMS11{s}';
end
save sim_data.mat simT simStraj simParams data

%%
simEx = cell(size(simStraj));
th = cell(size(simStraj));
h = waitbar(0,'looking for unspecific sticking events...');
for isp = 1:length(simStraj) %spot index
    waitbar(isp/length(simStraj),h,['looking for unspecific sticking events in trajectory ' ...
        num2str(isp) ' of ' num2str(length(simStraj)) '.']);
    straj = simStraj{isp};
    tmpRMS = simRMS11{isp}';

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
        states{2} = [S(updn==1) L(updn==1)];
        states{1} = [S(updn==-1) L(updn==-1)];

        % sort out by density of values below a certain threshold defined by total
        % percentage below ? remove and iterate until there is zero change
        % 

        P = 0.01;
        Q = 0.5;
        delta = 5;
        th{isp} = 0;
        seglims = [1 length(simStraj{1})];
        int_states = states;
        % pick out the ones contained in the segment
        for s = 1%:2
            int_states{s}(sum(int_states{s},2)<seglims(1),:) = [];
            int_states{s}(int_states{s}(:,1)>seglims(2),:) = [];
            if int_states{s}(1,1)<seglims(1)
                int_states{s}(1,2) = sum(int_states{s}(1,:))-seglims(1);
                int_states{s}(1,1) = seglims(1);
            end
            int_states{s}(end,2) = min(int_states{s}(end,2),seglims(2)-int_states{s}(end,1));

            go_on = 1;
            counter = 0;
            tmpEx = [];
            while go_on
                tmpI = setdiff(intersect(find(straj==s),seglims(1):seglims(2)),tmpEx);
                tmpV = sort(tmpRMS(tmpI));
                tmpP = max(1,floor(P*length(tmpV))-1);
                th{isp} = tmpV(tmpP);
                tmpD = zeros(length(tmpRMS),1);
                for n = 1:size(int_states{s},1)
                    for k = int_states{s}(n,1):sum(int_states{s}(n,:))
                        tmpF = intersect(k-delta:k+delta,int_states{s}(n,1):sum(int_states{s}(n,:)));
                        tmpD(k) = sum(tmpRMS(tmpF)<th{isp})/length(tmpF);
                    end
                end
                tmp = find(tmpD>Q & tmpRMS<th{isp});
                medopt = 0;
                if medopt
                    submed = find(tmpRMS < median(tmpRMS(tmpEx)));
                    submed = submed(submed>seglims(1) & submed<seglims(2));
                    tmp = unique([tmp; submed]);
                    %tmp = union(tmp,tmpI(tmpRMS(tmpI)<mean([th median(tmpRMS(tmpEx))])));
                end
                if all(ismember(tmp,tmpEx))
                    go_on = 0;
                else
                    tmpEx = unique([tmpEx; tmp]);
                end
                counter = counter + 1;
                disp(counter)
            end
            simEx{isp} = [simEx{isp} ; tmpEx];
        end
    end
end
close(h)
%% Display plots
figure
seglims = [1 length(simStraj{1})];
for isp = 1:length(simStraj)
    hold off
    tmpRMS = simRMS11{isp}';
    plot_twostateRMS(tmpRMS,simStraj{isp},seglims)
    hold on
    plot(seglims,th{isp}*[1 1],'k--')
    plot(simEx{isp},tmpRMS(simEx{isp}),'ko')
    set(gca,'Xlim',seglims,'Ylim',[0 2])
    pause
end
%%
RMSstate = cell(2,1);
for i = 1:length(simRMS11)
    tmpRMS = simRMS11{i}';
    for k = 1:2
        tmpIND = find(simStraj{i}==k);
        RMSstate{k} = [RMSstate{k};tmpRMS(tmpIND)];
    end
end
%%
RMSmax = 10;
for i = 1:2
    RMSstate{k} = RMSstate{k}(:,RMSstate{k}(1,:)<RMSmax);
end

%% Get density below 1% threshold for all bound intervals
P = 0.01;
globThresh = zeros(2,1);
for k = 1:2
    tmpV = sort(RMSstate{k});
    tmpP = max(1,floor(P*length(tmpV)));
    globThresh(k) = tmpV(tmpP);
end
densities = cell(length(simStraj),2);
areas = cell(length(simStraj),2);
h = waitbar(0,'');
for isp = 1:length(simStraj)
    waitbar(isp/length(simStraj),h,['getting density and area statistics ' ...
        num2str(isp) ' of ' num2str(length(simStraj)) '.']);
    straj = simStraj{isp};
    tmpRMS = simRMS11{isp}';
    tmpXY = simXY{isp};
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
        states{2} = [S(updn==1) L(updn==1)];
        states{1} = [S(updn==-1) L(updn==-1)];
    end
    for k = 1:2
        densities{isp,k} = zeros(size(states{k}));
        areas{isp,k} = zeros(size(states{k}));
        for i = 1:size(densities{isp,k},1)
            tmpI = states{k}(i,1):sum(states{k}(i,:))-1;
            tmpA = tmpRMS(tmpI)-globThresh(k);
            tmpB = tmpRMS(tmpI)<=globThresh(k);
            densities{isp,k}(i,1) = sum(tmpB)/length(tmpB);
            areas{isp,k}(i,1) = abs(sum(tmpB.*tmpA));
            Nmax = 100;
            if states{k}(i,2)>Nmax
                tmpA2 = zeros(1,length(tmpB)-Nmax+1);
                tmpB2 = zeros(1,length(tmpB)-Nmax+1);
                for j = 1:length(tmpB2)
                    tmpA2(j) = abs(sum(tmpA(j:j+Nmax-1).*tmpB(j:j+Nmax-1)));
                    tmpB2(j) = sum(tmpB(j:j+Nmax-1));
                end
                densities{isp,k}(i,2) = max(tmpB2)/Nmax;
                areas{isp,k}(i,2) = max(tmpA2);
            else
                densities{isp,k}(i,2) = densities{isp,k}(i,1);
                areas{isp,k}(i,2) = areas{isp,k}(i,1);
            end
        end
    end
end
close(h)

%% Check out distribution of densities

allD = cell(2,1);
allDmax = cell(2,1);
allA = cell(2,1);
allAmax = cell(2,1);
for isp = 1:size(densities,1)
    for k = 1:2
        if ~isempty(densities{isp,k})
            allD{k} = [allD{k};densities{isp,k}(:,1)];
            allDmax{k} = [allDmax{k};densities{isp,k}(:,2)];
            allA{k} = [allA{k};areas{isp,k}(:,1)];
            allAmax{k} = [allAmax{k};areas{isp,k}(:,2)];
        end
    end
end
%%

figure('Units','normalized','Position',[0 0 1 1])
for k = 1:2
    subplot(2,2,1+(k-1)*2)
    histogram(allD{k},0:0.01:1.1)
    title(['allD' num2str(k)],'Fontsize',16)
    subplot(2,2,2+(k-1)*2)
    histogram(allDmax{k},0:0.01:1.1)
    title(['allDmax' num2str(k)],'Fontsize',16)
    for i = 1:4
        subplot(2,2,i)
        ylim([0 100])
    end
end

%% cdfs
fractions = 0:0.01:1;
relcountsD = cell(2,1);
for k = 1:2
    relcountsD{k} = zeros(size(fractions));
    for i =1:length(relcountsD)
        relcountsD{k}(i) = sum(allDmax{k}>=fractions(i));
    end
    relcountsD{k} = relcountsD{k}/length(allDmax{k});
end


figure('Units','normalized','Position',[0 0 1 1])
for k = 1:2
    subplot(2,1,k)
    semilogx(fractions,relcountsD{k})
    title(['cum D, state ' num2str(k)], 'FontSize', 16)
    %set(gca,'Xdir','reverse')
end
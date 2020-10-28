%% Give percentages of wrongly assigned frames per trajectory
Fwrong = zeros(size(simStraj));
Nwrong = 0;
Nall_pre = 0;

for i = 1:length(simStraj)
    Fwrong(i) = sum(state_trajectories{i}~=simStraj{i})/length(simStraj{i});
    Nwrong = Nwrong + sum(state_trajectories{i}~=simStraj{i});
    Nall_pre = Nall_pre + length(simStraj{i});
end

figure
histogram(Fwrong)
title(['fraction of frames wrongly assigned, total: ' num2str(Nwrong/Nall_pre)], 'FontSize', 12)

%% Compare relative occupancies in trajectories pre- vs. post-evalutation

rOcc_pre = zeros(size(simStraj));
rOcc_post = zeros(size(simStraj));

rOccAll_pre = 0;
rOccAll_post = 0;
Nall_pre = 0;
Nall_post = 0;

for i = 1:length(simStraj)
    rOcc_pre(i) = sum(simStraj{i}==1)/length(simStraj{i});
    rOcc_post(i) = sum(state_trajectories{i}==1)/length(state_trajectories{i});
    rOccAll_pre = rOccAll_pre + sum(simStraj{i}==1);
    rOccAll_post = rOccAll_post + sum(state_trajectories{i}==1);
    Nall_pre = Nall_pre + length(simStraj{i});
    Nall_post = Nall_post + length(state_trajectories{i});
end

rOccAll_pre = rOccAll_pre/Nall_pre;
rOccAll_post = rOccAll_post/Nall_post;

figure
subplot(2,1,1)
histogram(rOcc_pre,20)
title(['rel. occupancy pre, total mean: ' num2str(rOccAll_pre)], 'FontSize', 12)
tmpLim1 = get(gca,'Xlim');
subplot(2,1,2)
histogram(rOcc_post,20)
title(['rel. occupancy post, total mean: ' num2str(rOccAll_post)], 'FontSize', 12)
tmpLim2 = get(gca,'Xlim');
for i = 1:2
    subplot(2,1,i)
    xlim([min(tmpLim1(1),tmpLim2(1)) max(tmpLim1(2),tmpLim2(2))])
end

%% Look for states that are missed completely

states_pre = cell(length(simStraj),2);
states_missed = cell(size(states_pre));

N_missed = zeros(size(states_missed));
N_post = zeros(size(states_missed));
N_pre = zeros(size(states_missed));

for i = 1:length(states_pre)
    states_pre(i,:) = getStates(simStraj{i})';
    for s = 1:2
        states_missed{i,s} = zeros(size(states_pre{i,s},1),1);
        for j = 1:size(states_missed{i,s},1)
            states_missed{i,s}(j) = all(state_trajectories{i}(states_pre{i,s}(j,1):(sum(states_pre{i,s}(j,:))-1))~=s);
        end
        N_missed(i,s) = sum(states_missed{i,s});
        N_post(i,s) = size(states_pre{i,s},1)-N_missed(i,s);
        N_pre(i,s) = size(states_pre{i,s},1);
    end
end

% figure
% for s = 1:2
%     subplot(1,2,s)
%     histogram(N_missed(:,s),-.5:max(N_missed(:,s))+.5)
%     title(['N_{missed} state ' num2str(s) ' --- N_{total}: ' num2str(sum(N_missed(:,s)))], 'FontSize', 12, 'Interpreter', 'tex')
% end

for s = 2:-1:1
    T_missed{s} = zeros(sum(N_missed(:,s)),1);
    T_post{s} = zeros(sum(N_post(:,s)),1);
    T_pre{s} = zeros(sum(N_pre(:,s)),1);
    counter_missed = 0;
    counter_post = 0;
    counter_pre = 0;
    for i = 1:length(states_pre)
        T_missed{s}(counter_missed+(1:N_missed(i,s))) = states_pre{i,s}(states_missed{i,s}==1,2);
        counter_missed = counter_missed+N_missed(i,s);
        T_post{s}(counter_post+(1:N_post(i,s))) = states_pre{i,s}(states_missed{i,s}==0,2);
        counter_post = counter_post+N_post(i,s);
        T_pre{s}(counter_pre+(1:N_pre(i,s))) = states_pre{i,s}(:,2);
        counter_pre = counter_pre+N_pre(i,s);
    end
end

% figure
% for s = 1:2
%     subplot(2,2,s)
%     h = histogram(T_missed{s});
%     title(['T_{missed} state ' num2str(s) ' --- N_{total}: ' num2str(sum(N_missed(:,s)))], 'FontSize', 12, 'Interpreter', 'tex')
%     subplot(2,2,s+2)f
%     histogram(T_post{s},h.BinEdges)
%     title(['T_{post} state ' num2str(s) ' --- N_{total}: ' num2str(sum(N_post(:,s)))], 'FontSize', 12, 'Interpreter', 'tex')
%     xlim(h.Parent.XLim)
% end

T50 = [0 0];
figure('Units', 'normalized', 'Position', [0 0 1 1]) 
for s = 1:2
    subplot(1,2,s)
    if ~isempty(T_missed{s})
        edges = 1:ceil(max(T_missed{s}));
        h = [histcounts(T_missed{s},edges)' histcounts(T_post{s},edges)'];
        Ntot = sum(h(:));
        div = sum(h,2);
        for i = 1:2
            h(:,i) = h(:,i)./div;
        end
        T50(s) = edges(find(h(:,1)>0.5,1,'last')+1);
        bar(edges(1:end-1)+0.5,h,'stacked')  
        title(['T_{missed} and T_{post} state ' num2str(s) ' --- N_{total}: ' num2str(Ntot)], 'FontSize', 12, 'Interpreter', 'tex')
    else
        title(['No missed events for state ' num2str(s)], 'FontSize', 12, 'Interpreter', 'tex')
    end
end

save missT.mat T_pre T_post T_missed T50
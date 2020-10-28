%index_discard = []; %THIS SHOULD ALREADY EXIST AND BE CORRECT
stickinTervals = cell(size(state_trajectories));
stickinDices = cell(size(state_trajectories));
h = waitbar(0,'');
for isp = setdiff(1:size(stickinTervals,1),index_discard)
    if ~isempty(segments{isp})
        tmpRMS = data{indicesHMM(isp,1)}{indicesHMM(isp,2),1}.vwcm.rms10(arxv{isp}.segments(1):arxv{isp}.segments(end))';
        tmpXY = arxv{isp}.XY;
        tmpS = state_trajectories{isp};
        [stickinTervals{isp}, stickinDices{isp}] = get_sticky(tmpRMS, tmpXY, tmpS, segments{isp}, segmInds{isp}, globThreshs, cutoffD, medSigmas, limitSEM, gThreshs2);
    end
    waitbar(isp/(size(stickinTervals,1)),h,[num2str(100*round(isp/(size(stickinTervals,1)),2)) ' percent done.'])
end
close(h)
% %% Plot result
% cf = figure('Units', 'normalized', 'Position', [0 0 1 1]);
% subplot(3,1,1)
% plot(tmpRMS)
% plot(stickinDices,tmpRMS(stickinDices)
% subplot(3,1,2)
% plot(state_trajectories{isp})
% ylim([0 3])
% EX = zeros(size(state_trajectories{isp}));
% EX(stickinDices) = 1;
% subplot(3,1,3)
% plot(EX);
% ylim([-1 2])
% pause
% subplot(3,1,1)
% %xlim([arxv{isp}.segments(1) arxv{isp}.segments(end)])
% xlim([1 length(state_trajectories{isp})])
% subplot(3,1,2)
% xlim([1 length(state_trajectories{isp})])
% subplot(3,1,3)
% xlim([1 length(state_trajectories{isp})])
% pause
% close(cf)
%% Save data
save sticky.mat stickinTervals stickinDices

%% GUI for inspection of state-assigned trajectories:

%

%mlmodel = model8_7;
%xyHMM = arxv8_7.xyHMM;
%inDisp = Arxv.stickinDices;
inData = setdiff(1:size(stickinTervals,1),index_discard);
inDisp = indicesHMM(inData,:);
m_start = 1;
i = 1;%find(inDisp(:,1)==m_start,1);
while isempty(state_trajectories{i}) && i<=length(state_trajectories)
    i = i+1;
end
ts = figure('Units','normalized','Position',[0 0 1 1]);
for p = 1:4
    subplot(4,1,p)
end
bBack = uicontrol('Style', 'pushbutton', 'String', 'Back', 'Units', 'normalized', 'Position', [0.025 0.8 0.05 0.04], 'Callback', 'if i > 1 i = i-1; end, while isempty(state_trajectories{inData(i)}) && i>1 i = i-1; end, uiresume', 'FontSize', 12);
bNext = uicontrol('Style', 'pushbutton', 'String', 'Next', 'Units', 'normalized','Position', [0.925 0.8 0.05 0.04], 'Callback', 'if inData(i) < length(state_trajectories) i = i+1; end, while isempty(state_trajectories{inData(i)}) && i<length(state_trajectories) i = i+1; end, uiresume', 'FontSize', 12);
loLim = uicontrol('Style', 'edit', 'Units', 'normalized', 'Position', [0.025 0.2 0.03 0.03]);
hiLim = uicontrol('Style', 'edit', 'Units', 'normalized', 'Position', [0.06 0.2 0.03 0.03]);
bSet = uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'String', 'Set Xlims', 'Position', [0.025 0.15 0.065 0.04], 'Callback', 'for p = 1:4 subplot(4,1,p), xlim([str2double(loLim.String) str2double(hiLim.String)]); end', 'FontSize', 12);
bReset = uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'String', 'Reset', 'Position', [0.025 0.1 0.065 0.04], 'Callback', 'for p = 1:4 subplot(4,1,p), xlim auto; end', 'FontSize', 12);
bDone = uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'String', 'Done', 'Position', [0.925 0.15 0.05 0.04], 'Callback', 'go_on = 0; uiresume', 'FontSize', 12);
onoff = {'off','on'};
bShow = uicontrol('Style', 'togglebutton', 'String', 'Show/Hide', 'Units', 'normalized', 'Position', [0.925 0.25 0.065 0.04], ... 
                'Callback', 'for p = 1:numel(rmp) if ~isempty(rmp{p}) rmp{p}.Visible = onoff{mod(bShow.Value,2)+1}; end; end;', 'FontSize', 12);

W = 11;
szRmv = 4;
color = {'g','m'};
go_on = 1;
rmp = cell(4,2);
while go_on
    tmpXY = arxv{inData(i)}.XY;
    tmpS = state_trajectories{inData(i)};
    tmpRMS = data{inDisp(i,1)}{inDisp(i,2),1}.vwcm.rms10(arxv{inData(i)}.segments(1):arxv{inData(i)}.segments(end));
    %tmpRMS = data{inDisp(i,1)}{inDisp(i,2),1}.vwcm.rms10(segments{inData(i)}(1):segments{inData(i)}(end));
    plot_twostate(tmpXY,tmpS,tmpRMS');
    subplot(4,1,1)
    hold on
    plot(double(tmpS)+.75, 'k')
    offset = segments{inData(i)}(1)-1;
    for intidx = 1:size(segments{inData(i)},1)
        plot(segments{inData(i)}(intidx,:)-offset,gThreshs2(2,segmInds{inData(i)}(intidx))*[1 1],'b--')
        for k = 1:2      
            plot(segments{inData(i)}(intidx,:)-offset,globThreshs(k,segmInds{inData(i)}(intidx))*[1 1],'k--')
        end
    end
    ylim([0 3])
    title(['Movie ' num2str(inDisp(i,1)) ', spot ' num2str(inDisp(i,2)), ', index ' num2str(i) '/' num2str(length(inData))],'FontSize',16)
    tmpMedXY = zeros(size(tmpXY));
    for k = 1:2
        tmpMedXY(k,:) = meanfilt1_trunc(tmpXY(k,:),W);
        subplot(4,1,k+2)
        hold on
        plot(tmpMedXY(k,:),'Color',[1 1 1]*.6)
    end
    subplot(4,1,2)
    hold on
    plot(sqrt(tmpMedXY(1,:).^2+tmpMedXY(2,:).^2),'Color',[1 1 1]*.6)
    for K = 1:2
        for intidx = 1:size(arxv{inData(i)}.segments,1)
            plot(segments{inData(i)}(intidx,:)-offset,3/sqrt(W)*mean(medSigmas{K}(segmInds{inData(i)}(intidx),:))*[1 1],'k--')
        end
    end
    for K = 1:2
        tmpF = find(tmpS==K);
        tmpI = intersect(tmpF,stickinDices{inData(i)});
        if ~isempty(tmpI)           
            subplot(4,1,1)
            rmp{1,K} = plot(tmpI,tmpRMS(tmpI),[color{K} 'o'],'MarkerSize',szRmv);
            subplot(4,1,2)
            rmp{2,K} = plot(tmpI,sqrt(tmpXY(1,tmpI).^2+tmpXY(2,tmpI).^2),[color{K} 'o'],'MarkerSize',szRmv);
            for p = 1:2
                subplot(4,1,p+2)
                hold on
                rmp{p+2,K} = plot(tmpI,tmpXY(p,tmpI),[color{K} 'o'],'MarkerSize',szRmv);
            end
        else
            rmp = cell(4,2);
        end
    end
    bShow.Value = 1;
    uiwait(gcf)
end
close(ts)
%}
display('Done.')

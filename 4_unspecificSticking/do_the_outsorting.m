cutoffD = [0.4 0.15];
limitSEM = [4 4];
%{
%RESET
stateFrames = stateFramesZ;
densities = densitiesZ;
noneAboveP2 = noneAboveP2Z;
RMSintSeg = RMSintSegZ;
removInds = removIndsZ;
[allD, allDmax] = get_allDmax(densities, index_discard);
%}
GO_ON = 1;
while GO_ON
    % identify suspects for unspecific sticking -> to be removed
    %Nsaved = length(allDmax{1}) + length(allDmax{2});
    for K = 1:2 % SELECT STATE
        go_on = 1;
        while go_on
            %distlBelow = cell(2,1);
            %distXYBelow = cell(2,1);
            %indsBelow = cell(2,1);
%             for i = 1:2
%                 distXYBelow{i} = zeros(2,0);
%                 indsBelow{i} = zeros(0,3);
%             end
            %cutoffL = [20 20];
            %removInds = cell(length(medI),2);
            display(['Total number of non-discarded states is: ' num2str(length(allD{K}))])
            if any(allDmax{K}>cutoffD(K))
                display(['number of states with density above ' num2str(cutoffD(K)) ...
                                ' is ' num2str(sum(allDmax{K}>cutoffD(K)))])
                    for isp = setdiff(1:size(densities,1),index_discard)
                        if ~isempty(densities{isp,K})
                            tmpI = find((stateFrames{isp,K}(:,2)>((K==2)*(Nmax+2*floor(W/2))) & densities{isp,K}(:,2)>cutoffD(K)) ...
                                        | (K==2 & (densities{isp,K}(:,2)>=0.5 | noneAboveP2{isp,K} == 1)));         
                            if ~isempty(tmpI)
                                %disp(length(tmpI));
                                densities{isp,K}(tmpI,:) = NaN;
                                removInds{K}{isp,1} = [removInds{K}{isp,1} ; tmpI];
                                for i = 1:length(tmpI)
                                    removInds{K}{isp,2} = [removInds{K}{isp,2} ; (stateFrames{isp,K}(tmpI(i),1):sum(stateFrames{isp,K}(tmpI(i),:))-1)'];
                                end
                            end
                            for i = 1:2
                                removInds{K}{isp,i} = unique(removInds{K}{isp,i});
                            end
                            stateFrames{isp,K}(removInds{K}{isp,1},:) = NaN;
                        end
                    end


                %% Remove from global cell array for segments and states
                [RMSintSeg, xySeg] = indRemover(RMSintSeg, xySeg, removInds);
    %             for j = 1:size(RMSintSeg,2)
    %                 for isp = 1:size(removInds{K},1)
    %                     RMSintSeg{K,j}(:,RMSintSeg{K,j}(3,:)==isp & ismember(RMSintSeg{K,j}(4,:),removInds{K}{isp,2})) = [];
    %                     xySeg{K,j}(:,xySeg{K,j}(3,:)==isp & ismember(xySeg{K,j}(4,:),removInds{K}{isp,2})) = [];
    %                 end 
    %             end

                %% Get density below 0.01% threshold for all bound intervals
                [globThreshs, gThreshs2] = get_globThreshs(RMSintSeg, 0.01, 0.025);
    %             P = 0.01;
    %             for i = 1:size(globThreshs,2)
    %                 for k = 1:2
    %                     tmpV = sort(RMSintSeg{k,i}(1,:));
    %                     tmpP = max(1,floor(P*length(tmpV)));
    %                     globThreshs(k,i) = tmpV(tmpP);
    %                 end
    %             end
                h = waitbar(0,'');
                for isp = setdiff(1:size(densities,1),index_discard)
                    waitbar(isp/length(state_trajectories),h,['getting density and area statistics ' ...
                        num2str(isp) ' of ' num2str(length(state_trajectories)) '.']);
                    if ~isempty(segments{isp}) 
                        tmpRMS = data{indicesHMM(isp,1)}{indicesHMM(isp,2),1}.vwcm.rms10(arxv{isp}.segments(1):arxv{isp}.segments(end))';
                        tmpXY = arxv{isp}.XY;
                        tmpSeg = 1;
                        for i = 1:size(densities{isp,K},1)
                            if ~isnan(stateFrames{isp,K}(i,1))
                                while stateFrames{isp,K}(i,1) > segments{isp}(tmpSeg,2)
                                    tmpSeg = tmpSeg + 1;
                                end
                                tmpA = 1;
                                tmpB = zeros(1,stateFrames{isp,K}(i,2));
                                tmpI = stateFrames{isp,K}(i,1):sum(stateFrames{isp,K}(i,:))-1;
                                tmpF = 0;
                                while sum(stateFrames{isp,K}(i,:)) > segments{isp}(tmpSeg,2)
                                    tmpA = tmpA && ...
                                        all(tmpRMS(tmpI((tmpF+1):(segments{isp}(tmpSeg,2)-stateFrames{isp,K}(i,1))))<gThreshs2(K,segmInds{isp}(tmpSeg)));
                                    tmpB((tmpF+1):(segments{isp}(tmpSeg,2)-stateFrames{isp,K}(i,1))) = ...
                                        tmpRMS(tmpI((tmpF+1):(segments{isp}(tmpSeg,2)-stateFrames{isp,K}(i,1))))<=globThreshs(K,segmInds{isp}(tmpSeg));
                                    tmpF = segments{isp}(tmpSeg,2)-stateFrames{isp,K}(i,1);
                                    tmpSeg = tmpSeg + 1;
                                end
                                if tmpF < stateFrames{isp,K}(i,2)
                                    tmpA = tmpA && all(tmpRMS(tmpI(tmpF+1:end))<gThreshs2(K,segmInds{isp}(tmpSeg)));
                                    tmpB(tmpF+1:end) = tmpRMS(tmpI(tmpF+1:end))<=globThreshs(K,segmInds{isp}(tmpSeg));
                                end
                                densities{isp,K}(i,1) = sum(tmpB)/length(tmpB);
                                noneAboveP2{isp,K}(i) = tmpA;
                                Nmax = 100;
                                Lmin = [Nmax, (Nmax+2*floor(W/2))];
                                if stateFrames{isp,K}(i,2)>Lmin(k)
                                    if K==2
                                        tmpB = tmpB(floor(W/2)+1:end-floor(W/2));
                                    end
                                    tmpB2 = zeros(1,length(tmpB)-Nmax+1);
                                    for j = 1:length(tmpB2)
                                        tmpB2(j) = sum(tmpB(j:j+Nmax-1));
                                    end
                                    densities{isp,K}(i,2) = max(tmpB2)/Nmax;
                                elseif K==1
                                    densities{isp,K}(i,2) = densities{isp,K}(i,1);
                                elseif K==2
                                    if densities{isp,K}(i,1)>=0.5
                                        densities{isp,K}(i,2) = densities{isp,K}(i,1);
                                    else
                                        densities{isp,K}(i,2) = NaN;
                                    end
                                end
                            else
                                densities{isp,K}(i,:) = NaN; % REDUNDANCY
                            end
                        end

                    end
                end
                close(h)

                %%
                [allD, allDmax] = get_allDmax(densities, index_discard);
        %         allD = cell(2,1);
        %         allDmax = cell(2,1);
        %         for k = 1:2
        %             for isp = setdiff(1:size(densities,1),index_discard)
        %                 if ~isempty(densities{isp,k})
        %                     allD{k} = [allD{k};densities{isp,k}(~isnan(densities{isp,k}(:,1)),1)];
        %                     allDmax{k} = [allDmax{k};densities{isp,k}(~isnan(densities{isp,k}(:,2)),2)];
        %                 end
        %             end
        %         end
            else
                display(['No states with density above ' num2str(cutoffD(K))])
                go_on = 0;
            end
            display(['Total number of non-discarded states is: ' num2str(length(allD{K}))])
        end
    end
    W = 11;
    suspects = cell(size(RMSintSeg));
    suspectSpots = [];
    suspectsRMS = cell(size(RMSintSeg));
    suspectsOveRel = cell(size(RMSintSeg));
    for i = numel(suspectsRMS):-1:1
        suspectsRMS{i} = zeros(3,0);
        suspectsOveRel{i} = zeros(3,0);
    end
    for K = 1:2
        for intidx = 1:size(RMSintSeg,2)
        %     suspects{intidx} = union(suspects{intidx},find(RMSintSeg{K,intidx}(1,:)<globThreshs(K,intidx) & ...
        %                                 (3/sqrt(W)*abs(xySeg{K,intidx}(1,:))>medSigmas{K}(intidx,1) | ...
        %                                 3/sqrt(W)*abs(xySeg{K,intidx}(2,:))>medSigmas{K}(intidx,2))));
            tmpI = find(RMSintSeg{K,intidx}(1,:)<globThreshs(K,intidx) & ...
                    sqrt((xySeg{K,intidx}(1,:)).^2 + (xySeg{K,intidx}(2,:)).^2)>limitSEM(K)/sqrt(W)*mean(medSigmas{K}(intidx,:)));
            if ~isempty(tmpI)
                suspects{K,intidx} = union(suspects{K,intidx},tmpI);
                suspectSpots = union(suspectSpots,unique(xySeg{K,intidx}(3,suspects{K,intidx})));
                suspectsRMS{K,intidx} = [suspectsRMS{K,intidx} RMSintSeg{K,intidx}([1 3 4],tmpI)];
                suspectsOveRel{K,intidx} = [suspectsOveRel{K,intidx} [sqrt((xySeg{K,intidx}(1,tmpI)).^2 + (xySeg{K,intidx}(2,tmpI)).^2).*sqrt(W)./mean(medSigmas{K}(intidx,:)); ...
                                                                    xySeg{K,intidx}(3:4,tmpI)]];
            end
        end
    end

    % Lengths of intervals fulfilling criteria
    suspectsL = cell(size(suspectsRMS));
    for k = 1:2
        foo = cat(2,suspectsRMS{k,:});
        foo(1,:) = 1:size(foo,2);
        foo = sortrows(foo', [2 3])';

        a = 1;
        i = 1;
        l = 1;
        L = zeros(1,size(foo,2));
        while i<=size(foo,2)-1
            if foo(2,i+1)==foo(2,i) && foo(3,i+1)==foo(3,i)+1
                l = l+1;
            else
                L(a:i) = l;
                a = i+1;
                l = 1;
            end
            i = i+1;
        end
        if a==i
            L(i) = 1;
        else
            L(a:i) = l;
        end

        [~,I] = sort(foo(1,:));
        L = L(I);

        offset = 0;
        for i = 1:size(suspectsL,2)
            suspectsL{k,i} = L((1:size(suspectsRMS{k,i},2))+offset);
            offset = offset + size(suspectsRMS{k,i},2);
        end
    end

    suspectInds = cell(size(removInds));
    for K = 1:2
        suspectInds{K} = cell(size(removInds{K}));
        for isp = reshape(suspectSpots,1,[]) 
            tmp = zeros(0,1);
            for intidx = 1:size(suspects,2)
                tmpI = sort(RMSintSeg{K,intidx}(4,suspects{K,intidx}(RMSintSeg{K,intidx}(3,suspects{K,intidx})==isp)));
                i = 1;
                while i <= numel(tmpI)
                    tmp = union(tmp,(tmpI(i)+(-2*floor(W/2):2*floor(W/2)))');
                    i = i+1;
                end
            end
            if ~isempty(tmp)
                suspectInds{K}{isp,2} = sort(unique(tmp));
            end
            tmp = [];
            for i = 1:size(stateFrames{isp,K},1)
                if ~isnan(stateFrames{isp,K}(i,1))
                    if ~isempty(intersect(stateFrames{isp,K}(i,1):sum(stateFrames{isp,K}(i,:)),suspectInds{K}{isp,2}))
                        tmp = [tmp; i];
                    end
                end
            end
            suspectInds{K}{isp,1} = tmp;
            %suspectFrames{i} = sort(RMSintSeg{K,intidx}(4,suspects(RMSintSeg{K,intidx}(3,suspects)==i)));
        end
    end

    %% Remove suspects from RMSintSeg, update density statistics:
    [RMSintSeg, xySeg, Nremoved] = indRemover(RMSintSeg, xySeg, suspectInds);
    display(['Total number of indices removed: ' num2str(Nremoved)])
    GO_ON = Nremoved>0;
    if GO_ON
        %%
        [globThreshs, gThreshs2] = get_globThreshs(RMSintSeg, 0.01, 0.025);
        %%
        h = waitbar(0,'');
        for isp = setdiff(1:size(densities,1),index_discard)
            for K = 1:2
                waitbar(isp/length(state_trajectories),h,['getting density and area statistics ' ...
                    num2str(isp) ' of ' num2str(length(state_trajectories)) '.']);
                if ~isempty(segments{isp})
                    tmpRMS = data{indicesHMM(isp,1)}{indicesHMM(isp,2),1}.vwcm.rms10(arxv{isp}.segments(1):arxv{isp}.segments(end))';
                    Lmin = [Nmax, (Nmax+2*floor(W/2))];
                    [densities{isp,K}, noneAboveP2{isp,K} ] = update_density(densities{isp,K}, noneAboveP2{isp,K}, tmpRMS, globThreshs(K,:), segments{isp}, ...
                                                                            segmInds{isp}, stateFrames{isp,K}, Lmin(K), gThreshs2(K,:));
                end
            end
        end
        close(h)
        %%
        [allD, allDmax] = get_allDmax(densities, index_discard);
    end
end
display('Established final estimate for global thresholds.')

%
%execute_sticky
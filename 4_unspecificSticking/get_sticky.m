function [ intervals, indices ] = get_sticky( RMS, XY, straj, segments, segmInds, threshs, cutoffD, sigmas, multSEM, threshs2 )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
    W = 11;
    Nmax = 100;
    Lmin = [Nmax, (Nmax+2*floor(W/2))];
    %cutoffD = [.25 .5];
    if nargin < 10
        threshs2 = zeros(size(threshs));
    end
    RMS = reshape(RMS,1,[]);
    
    SEM3 = cell(2,1);  
    for i = 1:2
        XY(i,:) = meanfilt1_trunc(XY(i,:),W);
        SEM3{i} = zeros(size(sigmas{i},1),1);
        for j = 1:size(sigmas{i},1)
            SEM3{i}(j) = multSEM(i)*mean(sigmas{i}(j,:))./sqrt(W);
        end
    end
    rmXY = sqrt(XY(1,:).^2+XY(2,:).^2);
    
    exclude = zeros(size(straj));
    
    %find states
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
            states{k}(sum(states{k},2)>segments(end),:) = [];
        end
    end
    %disp(states{1})
    %disp(states{2})
    
    %exclude by density criterion
    %
    for k = 2%1:2
        tmpSeg = 1;
        for i = 1:size(states{k},1)
            while states{k}(i,1) > segments(tmpSeg,2)
                tmpSeg = tmpSeg + 1;
            end
            tmpA = states{k}(i,2)>1;
            tmpB = zeros(1,states{k}(i,2));
            tmpI = states{k}(i,1):sum(states{k}(i,:))-1;
            tmpF = 0;
            while sum(states{k}(i,:)) > segments(tmpSeg,2)
                tmpA = tmpA && ...
                    all(RMS(tmpI((tmpF+1):(segments(tmpSeg,2)-states{k}(i,1))))<threshs2(k,segmInds(tmpSeg)));
                tmpB((tmpF+1):(segments(tmpSeg,2)-states{k}(i,1))) = ...
                    RMS(tmpI((tmpF+1):(segments(tmpSeg,2)-states{k}(i,1))))<=threshs(k,segmInds(tmpSeg));
                tmpF = segments(tmpSeg,2)-states{k}(i,1);
                tmpSeg = tmpSeg + 1;
            end
            if tmpF < states{k}(i,2)
                tmpA = tmpA && all(RMS(tmpI(tmpF+1:end))<threshs2(k,segmInds(tmpSeg)));
                tmpB(tmpF+1:end) = RMS(tmpI(tmpF+1:end))<=threshs(k,segmInds(tmpSeg));
            end
            if tmpA==1 && k==2
                exclude(states{k}(i,1)-1:sum(states{k}(i,:))) = exclude(states{k}(i,1)-1:sum(states{k}(i,:))) + 1;
            elseif any(tmpB==1)
                if states{k}(i,2)>Lmin(k)
                    if k==2
                        tmpB = tmpB(floor(W/2)+1:end-floor(W/2));
                    end
                    tmpB2 = zeros(1,length(tmpB)-Nmax);
                    for j = 1:length(tmpB2)
                        tmpB2(j) = sum(tmpB(j:j+Nmax));
                    end
                    tmpB2 = tmpB2./Nmax;
                    tmpE = zeros(size(tmpB));
                    tmpE(Nmax/2+1:end-Nmax/2) = tmpB2>cutoffD(k);
                    tmpE(1:Nmax/2) = tmpE(Nmax/2+1);
                    tmpE(end-Nmax/2+1:end) = tmpE(end-Nmax/2);
                    if k==2
                        tmpE = [tmpE(1)*ones(1,floor(W/2)) tmpE tmpE(end)*ones(1,floor(W/2))];
                    end
                    exclude(states{k}(i,1):sum(states{k}(i,:))-1) = exclude(states{k}(i,1):sum(states{k}(i,:))-1) + tmpE;
                    if tmpE(1)==1
                        exclude(states{k}(i,1)-1) = exclude(states{k}(i,1)-1) + 1;
                    end
                    if tmpE(end)==1
                        exclude(sum(states{k}(i,:))) = exclude(sum(states{k}(i,:))) + 1;
                    end
                elseif k==2 && states{k}(i,2)>1 && (sum(tmpB)/length(tmpB) >= 0.5)
                    exclude(states{k}(i,1)-1:sum(states{k}(i,:))) = exclude(states{k}(i,1)-1:sum(states{k}(i,:))) + 1;
                end
            end
        end         
    end
    
    %}
    
    %exclude by off-center criterion and if RMS in state 2 is lower than threshold for state 1
    tmpE = zeros(size(straj));
    for k = 1:2
        tmpF = find(straj==k);
        for i = 1:size(segments,1)
            tmpE(tmpF(ismember(tmpF,segments(i,1):segments(i,2)) & ...
                RMS(tmpF)<threshs(k,segmInds(i)) & ...
                rmXY(tmpF)>SEM3{k}(segmInds(i)))) = 1;
            if k==2
                tmpE(tmpF(ismember(tmpF,segments(i,1):segments(i,2)) & ...
                RMS(tmpF)<threshs(1,segmInds(i)))) = 1;
            end
        end
    end
    for i = find(tmpE==1)
        tmpE(max(1,i-floor(W/2)):min(length(tmpE),i+floor(W/2))) = 1;
    end
    exclude = exclude + tmpE;
    %}
    
    %transform to intervals and indices
    if any(exclude>0)
        exclude = exclude>0;
        steps = find(diff(exclude)~=0) + 1;
        % Start frames and lengths of intervals
        S = reshape(steps(1:end-1),length(steps)-1,1);
        L = reshape(steps(2:end)-steps(1:end-1),length(steps)-1,1);
        % Assign type of interval
        updn = zeros(size(S));
        for i = 1:length(updn)
            updn(i) = sign(exclude(steps(i))-exclude(steps(i)-1));
        end
        % write excluded intervals
        intervals = [S(updn==1) L(updn==1)];
        intervals(:,2) = sum(intervals,2)-1;
        indices = find(exclude==1);
    else
        intervals = zeros(0,2);
        indices = [];
    end
end


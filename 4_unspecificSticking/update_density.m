function [ density, noneAboveP2 ] = update_density( density, noneAboveP2, RMS, threshs, segments, segmInds, stateFrames, Lmin, threshs2 )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    tmpSeg = 1;
    for i = 1:size(density,1)
        if ~isnan(stateFrames(i,1))
            while stateFrames(i,1) > segments(tmpSeg,2)
                tmpSeg = tmpSeg + 1;
            end
            tmpA = 1;
            tmpB = zeros(1,stateFrames(i,2));
            tmpI = stateFrames(i,1):sum(stateFrames(i,:))-1;
            tmpF = 0;
            while sum(stateFrames(i,:)) > segments(tmpSeg,2)
                tmpA = tmpA && all(RMS(tmpI((tmpF+1):(segments(tmpSeg,2)-stateFrames(i,1))))<threshs2(segmInds(tmpSeg)));
                tmpB((tmpF+1):(segments(tmpSeg,2)-stateFrames(i,1))) = ...
                    RMS(tmpI((tmpF+1):(segments(tmpSeg,2)-stateFrames(i,1))))<=threshs(segmInds(tmpSeg));
                tmpF = segments(tmpSeg,2)-stateFrames(i,1);
                tmpSeg = tmpSeg + 1;
            end
            if tmpF < stateFrames(i,2)
                tmpA = tmpA && all(RMS(tmpI(tmpF+1:end))<threshs2(segmInds(tmpSeg)));
                tmpB(tmpF+1:end) = RMS(tmpI(tmpF+1:end))<=threshs(segmInds(tmpSeg));
            end
            density(i,1) = sum(tmpB)/length(tmpB);
            noneAboveP2(i) = tmpA;
            Nmax = 100;
            if stateFrames(i,2)>Lmin
                if Lmin>Nmax
                    W = Lmin-Nmax;
                    tmpB = tmpB(floor(W/2)+1:end-floor(W/2));
                end
                tmpB2 = zeros(1,length(tmpB)-Nmax+1);
                for j = 1:length(tmpB2)
                    tmpB2(j) = sum(tmpB(j:j+Nmax-1));
                end
                density(i,2) = max(tmpB2)/Nmax;
            elseif ~isnan(density(i,2));
                density(i,2) = density(i,1);
            end
        else
            density(i,:) = NaN; % REDUNDANCY
        end
    end
end


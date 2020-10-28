function [ xyOut, meansCell, nCell ] = correctXYensemble( xyIn, indices, intervals )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

M = max(indices(:,1));
xyOut = cell(size(xyIn));
meansCell = cell(M,1);
nCell = cell(M,1);
counter = 0;
for m = 1:M
    tmpIntervals = intervals(indices(:,1)==m,:);
    zeroFrame = min(tmpIntervals(:,1)) - 1;
    tmpMean = zeros(2,max(tmpIntervals(:,2))); 
    tmpN = zeros(1,max(tmpIntervals(:,2)));
    tmpL = size(tmpIntervals,1);
    if tmpL>=3
        for i = 1:tmpL
            counter = counter + 1;
            tmpI = tmpIntervals(i,:);
            tmpMean(:,tmpI(1):tmpI(2)) = tmpMean(:,tmpI(1):tmpI(2)) + xyIn{counter};
            tmpN(tmpI(1):tmpI(2)) = tmpN(tmpI(1):tmpI(2)) + 1;
        end
        tmpMean = tmpMean(:,zeroFrame+1:end);
        tmpN = tmpN(zeroFrame+1:end);
        tmpMean = tmpMean./[tmpN;tmpN];
        counter = counter - tmpL;    
        for i = 1:tmpL
            counter = counter + 1;
            tmpI = tmpIntervals(i,:);
            xyOut{counter} = xyIn{counter} - tmpMean(:,(tmpI(1):tmpI(2))-zeroFrame);
        end
    else
        xyOut(counter+(1:tmpL)) = xyIn(counter+(1:tmpL));
        counter = counter+tmpL;
    end
    meansCell{m} = tmpMean;
    nCell{m} = tmpN;
end


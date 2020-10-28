function [RMSintSeg, xySeg, Nremoved] = indRemover(RMSintSeg, xySeg, removInds)
    Nremoved = 0;
    for K = 1:2
        for j = 1:size(RMSintSeg,2)
            tmp = size(RMSintSeg{K,j},2);
            for isp = 1:size(removInds{K},1)
                RMSintSeg{K,j}(:,RMSintSeg{K,j}(3,:)==isp & ismember(RMSintSeg{K,j}(4,:),removInds{K}{isp,2})) = [];
                xySeg{K,j}(:,xySeg{K,j}(3,:)==isp & ismember(xySeg{K,j}(4,:),removInds{K}{isp,2})) = [];
            end
            Nremoved = Nremoved + tmp - size(RMSintSeg{K,j},2);
        end
    end
end
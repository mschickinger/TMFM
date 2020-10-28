 function [ globThreshs, gThreshs2 ] = get_globThreshs( RMSintSeg, P, P2 )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if nargin==2
    P2 = 0.1;
end
    globThreshs = zeros(size(RMSintSeg));
    gThreshs2 = globThreshs;
    for i = 1:size(RMSintSeg,2)
        for k = 1:2
            if ~isempty(RMSintSeg{k,i})
                tmpV = sort(RMSintSeg{k,i}(1,:));
                tmpP = max(1,floor(P*length(tmpV)));
                globThreshs(k,i) = tmpV(tmpP);
                tmpP2 = max(1,floor(P2*length(tmpV)));
                gThreshs2(k,i) = tmpV(tmpP2);
            end
        end
    end
end


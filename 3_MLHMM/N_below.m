function [ output ] = N_below (traces,threshs )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N = zeros(1,numel(threshs));
N_all = 0;
threshs = sort(threshs);

if iscell(traces)
    for i = 1:length(traces)
        N_all = N_all + nnz(traces{i});
        tmp_trace = fliplr(reshape(nonzeros(traces{i}),1,[]));
        for j = 1:numel(threshs)
            tmpF = find(tmp_trace > threshs(j),1);
            if ~isempty(tmpF)
                N(j) = N(j) + tmpF-1;
            else
                N(j) = N(j) + numel(tmp_trace);
            end
        end
    end
else
    fprintf('Warning: Input _traces_ is not a cell array.\nTreating as numeric array with dimensions (N_traces)x(N_frames).\n')
    for i = 1:size(traces,1)
        N_all = N_all + nnz(traces(i,:));
        tmp_trace = fliplr(reshape(nonzeros(traces(i,:)),1,[]));
        for j = 1:numel(threshs)
            tmpF = find(tmp_trace > threshs(j),1);
            if ~isempty(tmpF)
                N(j) = N(j) + tmpF-1;
            else
                N(j) = N(j) + numel(tmp_trace);
            end
        end
    end
end
    
output.N = N;
output.N_all = N_all;
output.threshs = threshs;
end


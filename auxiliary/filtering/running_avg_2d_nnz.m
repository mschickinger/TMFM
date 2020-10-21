function [ means, varargout ] = running_avg_2d_nnz( trace, wSize )
% running_avg_2d_nnz: Returns local average positions for 2dim raw data 
% (eg estimated positions) ? only includes non-zero values in average

delta=floor(wSize/2);
varargout{1} = delta;
means=zeros(length(trace),2);
for i=1:delta
    tmp = trace(1:i+delta,:);
    tmp1 = tmp(tmp(:,1)>0,1);
    tmp2 = tmp(tmp(:,2)>0,2);
    means(i,:) = [mean(tmp1) mean(tmp2)];
end
for i=delta+1:length(trace)-delta
    tmp = trace(i-delta:i+delta,:);
    tmp1 = tmp(tmp(:,1)>0,1);
    tmp2 = tmp(tmp(:,2)>0,2);
    means(i,:) = [mean(tmp1) mean(tmp2)];
end
for i=length(trace)-delta:length(trace)
    tmp = trace(i-delta:end,:);
    tmp1 = tmp(tmp(:,1)>0,1);
    tmp2 = tmp(tmp(:,2)>0,2);
    means(i,:) = [mean(tmp1) mean(tmp2)];
end
end


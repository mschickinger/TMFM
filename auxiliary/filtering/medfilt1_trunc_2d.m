function [ output ] = medfilt1_trunc_2d( input, order )
%MEDFILT1_TRUNC median filter with endpoint correction
%   Corrects for "zero assumption" in medfilt1.
%   Truncates windows for median at start and end of the trace
sizeIn = size(input);
input = reshape(input,[length(input),2]);
output = [medfilt1(input(:,1),order) , medfilt1(input(:,2),order)];
for d = 1:2
    for i = 1:ceil((order-1)/2)
        output(i,d) = median(input(1:i+ceil((order-1)/2),d));
    end
    for i = length(input)-ceil((order-1)/2)+1:length(input)
        output(i,d) = median(input(i-ceil((order-1)/2):end,d));
    end
end
output = reshape(output, sizeIn);
end
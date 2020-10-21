function [ output ] = medfilt1_trunc( input, order )
%MEDFILT1_TRUNC median filter with endpoint correction
%   Corrects for "zero assumption" in medfilt1.
%   Truncates windows for median at start and end of the trace

output = medfilt1(input, order);
for i = 1:ceil((order-1)/2)
    output(i) = median(input(1:i+ceil((order-1)/2)));
end
for i = length(input)-ceil((order-1)/2)+1:length(input)
    output(i) = median(input(i-ceil((order-1)/2):end));
end
end
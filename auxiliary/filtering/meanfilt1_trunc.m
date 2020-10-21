function [ output ] = meanfilt1_trunc( input, order )
%MEANFILT1_TRUNC mean filter with endpoint correction.
%   Truncates windows for filtering at start and end of the trace

w = ceil((order-1)/2);

output = zeros(size(input));
for i = 1:w
    output(i) = mean(input(1:i+w));
end
for i = w+1:length(input)-w
    output(i) = mean(input(i-w:i+w));
end
for i = length(input)-w+1:length(input)
    output(i) = mean(input(i-ceil((order-1)/2):end));
end
end
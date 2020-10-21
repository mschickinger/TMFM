function [ logP ] = log_erlangpdf(x,lambda,n)
% n: number of events
% x: time for n events (or sum of n lifetimes)
% lambda: mean probability of occurence

logP = n.*log(lambda) + (n-1).*log(x) - sum(log(1:(n-1))) - lambda.*x;

end


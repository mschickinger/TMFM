function [ P ] = erlangcdf(x,lambda,n)
% n: number of events
% x: time for n events (or sum of n lifetimes)
% lambda: mean probability of occurence

P = gammainc(lambda.*x,n);

end


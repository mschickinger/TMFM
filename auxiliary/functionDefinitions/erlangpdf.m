function [ p ] = erlangpdf(x,lambda,n)
% n: number of events
% x: time for n events (or sum of n lifetimes)
% lambda: mean probability of occurence

p = (lambda.^n).*(x.^(n-1))./gamma(n).*exp(-lambda.*x);

end


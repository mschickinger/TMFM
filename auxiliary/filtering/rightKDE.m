function [ n, p, x_points] = rightKDE( x_data, h, n_points, varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    q = inputParser;
    addRequired(q, 'x_data')
    addRequired(q, 'h')
    addRequired(q, 'n_points')
    addParameter(q, 'kernel', 'box', @ischar)
    parse(q,x_data,h,n_points,varargin{:})
    
    x_data = q.Results.x_data;
    h = q.Results.h;
    n_points = q.Results.n_points;
    
    n = zeros(1,n_points);
    
    if strcmp(q.Results.kernel,'box')
        x_points = linspace(min(x_data)+h/2,max(x_data),n_points);
        delta = h/2;
        for i= 1:length(x_points) % loop over data
            n(i) = sum((x_data>=x_points(i)-delta) & (x_data<x_points(i)+delta))/h;
        end
    elseif strcmp(q.Results.kernel,'gauss')
        %display('Gauss, baby, gauss!')
        x_points = linspace(min(x_data)+3*h,max(x_data),n_points);
        sigma = h;
        for i = 1:length(x_data)
            n = n + normpdf(x_points,x_data(i),sigma);
        end
    end    
    p = n/numel(x_data);
end


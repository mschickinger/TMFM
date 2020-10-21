function [ ] = errBar( X, Y, E, style, color, varargin )
% Scatter / Line plot with nicer errorbars

p = inputParser;
addRequired(p, 'X')
addRequired(p, 'Y')
addRequired(p, 'E')
addOptional(p, 'style', 'o')
addOptional(p, 'color', [0 0 0])
addParameter(p,'LineWidth', 1)
addParameter(p,'MarkerSize', 8)
addParameter(p,'CapWidth',

parse(p, X, Y, E, varargin{:})

X = p.Results.X;
Y = p.Results.Y;
E = p.Results.E;

if any(p.Results.style == '-')
    linestyle = '-';
    if strfind(bla,'--')
        linestyle = '--';
    end
else
    linestyle = 0;
end

if any(p.Results.style == 'x')
    marker = 'x';
elseif any(p.Results.style == 'o')
    marker = 'o';
end

% prepare arrays for errobars

xerr = zeros(3*numel(X)-1,1);
yerr = zeros(size(xerr));

xerr(1:3:end) = X;
xerr(2:3:end) = X;
xerr(3:3:end) = X(2:end);
yerr(1:3:end) = Y+E;
yerr(2:3:end) = Y-E;
yerr(3:3:end) = NaN;

plot(xerr,yerr,'-','Color',p.Results.color)

end


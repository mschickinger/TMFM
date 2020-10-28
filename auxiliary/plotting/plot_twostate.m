function plot_twostate(XY, S, varargin)

% parse input
p = inputParser;

addRequired(p,'XY')
addRequired(p,'S')
addOptional(p,'RMS',[])
addOptional(p,'R',[])
addParameter(p,'wSize',11)

parse(p, XY, S, varargin{:})

XY = p.Results.XY;
X = XY(1,:);
Y = XY(2,:);

S = p.Results.S;

if isempty(p.Results.RMS)
    RMS = RMSfilt2d(XY',p.Results.wSize)';
else
    RMS = p.Results.RMS;
end

if isempty(p.Results.R)
    R = sqrt(X.^2+Y.^2);
else
    R = p.Results.R;
end

T{2} = 1:length(X);
T{1} = T{2};

transitions = find(diff(S)~=0)+1;
interS = S(transitions);
interT = transitions - 0.5;
interRMS = (RMS(transitions)+RMS(transitions-1))./2;

[mT{2}, IDX] = sort([T{2} interT]);
mT{1} = mT{2};
mS = [S interS];
mS = mS(IDX);
RMS = [RMS interRMS];
RMS = RMS(IDX);

s = mS(1);
mT{mod(s,2)+1}(1) = NaN;
i = 1;
while i<=length(mT{1})-1
    d = find(mS(i+1:end)~=s,1);
    if ~isempty(d)
        mT{mod(s,2)+1}(i+1:i+d-1) = NaN;
        i = i + d;
        s = mod(s,2) + 1;
    else
        mT{mod(s,2)+1}(i+1:end) = NaN;
        i = length(mT{1});
    end
end

T{1}(S==2) = NaN;
T{2}(S==1) = NaN;

% parameters for plotting
c = {[204 0 0]/255,[0 102 153]/255};
%c = {[0.85 0.325 0.098],[0 0.447 0.741]};
%c = {[204 153 0]/255,[0 102 153]/255};
%c = {'red','blue'};
mSize = 5;

subplot(4,1,1)
% RMS fluctuation
hold off
for s = 1:2
    plot(mT{s},RMS,'-','Color',c{s})
    hold on
end

subplot(4,1,2)
% radius
hold off
for s = 1:2
    plot(T{s},R,'.','Color',c{s},'MarkerSize', mSize)
    hold on
end

subplot(4,1,3)
% x-displacement
hold off
for s = 1:2
    plot(T{s},X,'.','Color',c{s},'MarkerSize', mSize)
    hold on
end

subplot(4,1,4)
% y-displacement
hold off
for s = 1:2
    plot(T{s},Y,'.','Color',c{s},'MarkerSize', mSize)
    hold on
end

end
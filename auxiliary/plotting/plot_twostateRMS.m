function plot_twostateRMS(mR, S, segments)

mR = mR(segments(1):segments(end))';
T{2} = 1:length(mR);
T{1} = T{2};

transitions = find(diff(S)~=0)+1;
interS = S(transitions);
interT = transitions - 0.5;
interR = (mR(transitions)+mR(transitions-1))./2;

[mT{2}, IDX] = sort([T{2} interT]);
mT{1} = mT{2};
mS = [S interS];
mS = mS(IDX);
mR = [mR interR];
mR = mR(IDX);

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

% RMS trace
hold off
for s = 1:2
    plot(mT{s},mR,'-','Color',c{s})
    hold on
end

end
function [ lifetimes, Nremoved ] = lts_strict_cutoff( state_trajectories, Fmin, varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

p = inputParser;

addRequired(p, 'state_trajectories')
addRequired(p, 'Fmin')
addOptional(p, 'ex_int', 0)

parse(p, state_trajectories, Fmin, varargin{:})

straj = p.Results.state_trajectories;
Fmin = p.Results.Fmin;
if iscell(p.Results.ex_int)
    ex_int = p.Results.ex_int;
else
    ex_int = cell(size(straj));
end
%tpf = p.Results.tpf;

lifetimes{2} = zeros(0,1);
lifetimes{1} = lifetimes{2};
Nremoved = [0 0];

for i = 1:length(straj)
    tmpS = getStates(straj{i},ex_int{i});
    for s = 1:2
        if ~isempty(tmpS{s})
            tmpS{s}(tmpS{s}(:,2)>=Fmin(s),:) = [];
            Nremoved(s) = Nremoved(s) + size(tmpS{s},1);
        end
    end
    tmpS = sortrows(vertcat(tmpS{:}));
    for j = 1:size(tmpS,1)
        straj{i}(tmpS(j,1):sum(tmpS(j,:))-1) = straj{i}(tmpS(j,1)-1);
    end
    tmpS = getStates(straj{i},ex_int{i});
    for s = 1:2
        if ~isempty(tmpS{s})
            lifetimes{s} = [lifetimes{s};tmpS{s}(:,2)];
        end
    end
end    
    
return


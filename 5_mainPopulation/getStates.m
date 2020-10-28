function [ states ] = getStates( straj, varargin )

p = inputParser;

addRequired(p, 'straj')
addOptional(p, 'ex_int', [])

parse(p,straj,varargin{:})
straj = p.Results.straj;
ex_int = p.Results.ex_int;
%
    steps = find(diff(straj)~=0) + 1;
    if ~isempty(steps)
        % Start frames and lengths of states
        S = reshape(steps(1:end-1),length(steps)-1,1);
        L = reshape(steps(2:end)-steps(1:end-1),length(steps)-1,1);
        % Assign type of states
        updn = zeros(size(S));
        for j = 1:length(updn)
            updn(j) = sign(straj(steps(j))-straj(steps(j)-1));
        end
        % divide in hi and lo states
        states{2} = [S(updn==1) L(updn==1)];
        states{1} = [S(updn==-1) L(updn==-1)];
        if ~isempty(ex_int)
            for k = 1:2
                states{k} = remove_excluded(states{k}, ex_int);
            end
        end
    else
        states = cell(2,1);
    end
    
    function [SL] = remove_excluded(SL, ex_int)
        tmp_ex = [];
        for i = 1:size(ex_int,1)
            tmp_ex = [tmp_ex ex_int(i,1):ex_int(i,2)];
        end
        i = 1;
        while i <= size(SL,1)
            if any(ismember(SL(i,1):(sum(SL(i,:))-1),tmp_ex))
                SL(i,:) = [];
            else
                i = i+1;
            end
        end
    end

end


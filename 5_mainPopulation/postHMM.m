function output = postHMM(INPUT, varargin)

%   'INPUT' is a struct with fields:
%       'indices' - (array  of movie and spot indices)
%       'ranges' - Frame ranges evaluated in the corresponding spot (start of first segment until end of last segment) 
%       'StateTrajectories' - Vector of 1s and 2s indicating the state assigned to corresponding frame
%       'XY' - XY-trajectories that were used during HMM evaluation
%       'medI' - median-filtered intensity traces of the particles.
%       'ex_int' - intervals to be excluded in further analysis.

%   'output' is a struct with fields:
%   'hop' - 
%   'scatterStats' - 
%   'allStats' - 
%   'indices' - 
%   'ranges' - 

p = inputParser;

addRequired(p,'INPUT');
addOptional(p,'sample_ident',[]);
addOptional(p,'tpf',[]);

parse(p,INPUT,varargin{:});

INPUT = p.Results.INPUT;

% Create hop structure
if ~isempty(p.Results.sample_ident)
    sample_ident = p.Results.sample_ident;
else
    sample_ident = inputdlg({'Date:', 'Sample:', 'Number of movies:'}, 'Identify');
end

hop.date = sample_ident{1};
hop.sample = sample_ident{2};
N_movies = str2double(sample_ident{3});

% Times per frame
if ~isempty(p.Results.tpf)
    hop.tpf = p.Results.tpf;
else
    input_lines = cell(N_movies,1);
    def_ans = cell(N_movies,1);
    for m = 1:N_movies
        input_lines{m} = ['Enter tpf for movie # ' num2str(m)];
        def_ans{m} = '50';
    end
    hop.tpf = str2double(inputdlg(input_lines, 'Times per frame', 1, def_ans));
end

% hop.results, scatterStats and allStats
hop.results = cell(N_movies,1);
scatterStats = zeros(size(INPUT.indices,1),6); % mTb mTu sDb sDu Nb Nu;
counter = 1;
N_hi = 0;
N_lo = 0;
for m = 1:N_movies
    tmp_spotnums = INPUT.indices(INPUT.indices(:,1)==m,2);
    hop.results{m} = cell(length(tmp_spotnums),1);
    for s = 1:length(tmp_spotnums)
        start_offset = INPUT.ranges(counter,1) - 1;
        hop.results{m}{s}.spotnum = tmp_spotnums(s);
        hop.results{m}{s}.state_trajectory = INPUT.state_trajectories{counter};
        hop.results{m}{s}.ex_int = INPUT.ex_int{counter};
        for k = 2:-1:1
            hop.results{m}{s}.relOcc(k) = relative_occupancy(INPUT.state_trajectories{counter}, k, INPUT.ex_int{counter}, start_offset);
        end
        [tmp_hi, tmp_lo] = get_hilo(hop.results{m}{s}.state_trajectory, start_offset, INPUT.ex_int{counter});
        if ~isempty(tmp_hi)    
            hop.results{m}{s}.hi = tmp_hi;
            scatterStats(counter,1) = mean(tmp_hi(:,2));
            scatterStats(counter,3) = std(tmp_hi(:,2));
        else
            hop.results{m}{s}.hi = [];
        end    
        if ~isempty(tmp_lo)
            hop.results{m}{s}.lo = tmp_lo;
            scatterStats(counter,2) = mean(tmp_lo(:,2));  
            scatterStats(counter,4) = std(tmp_lo(:,2));
        else
            hop.results{m}{s}.lo = [];
        end
        % conversion from frames to seconds:
        scatterStats(counter,1:4) = 2*hop.tpf(m)/1000*scatterStats(counter,1:4);
        % Number of states (hi, lo)
        scatterStats(counter,5) = size(tmp_hi,1);
        N_hi = N_hi + size(tmp_hi,1);
        scatterStats(counter,6) = size(tmp_lo,1);
        N_lo = N_lo + size(tmp_lo,1);
        counter = counter + 1;
    end
end

% remove all spots without transitions
% tmp_remove = find(scatterStats(:,1) == 0);
% tmp_keep = find(scatterStats(:,1) ~= 0);
% scatterStats(tmp_remove,:) = [];
% INPUT.XY(tmp_remove) = [];
% INPUT.medI(tmp_remove) = [];
hop.indices = INPUT.indices;%(tmp_keep,:);
hop.ranges = INPUT.ranges;%(tmp_keep,:);

% allspotStats
allStats.hi = zeros(N_hi,10);
allStats.lo = zeros(N_lo,10);
counterHi = 0;
counterLo = 0;
counter = 1;
for m = 1:N_movies
    for s = 1:size(hop.results{m},1)
        start_offset = INPUT.ranges(counter,1) - 1;
        tmpNhi = size(hop.results{m}{s}.hi,1);
        if tmpNhi>0
            %first column: movie index
            allStats.hi(counterHi+(1:tmpNhi),1) = m;
            %second column: spot index
            allStats.hi(counterHi+(1:tmpNhi),2) = hop.results{m}{s}.spotnum;
            %third/fourth column: start and duration (frames)
            allStats.hi(counterHi+(1:tmpNhi),3:4) = hop.results{m}{s}.hi;
            %fifth column: duration (seconds)
            allStats.hi(counterHi+(1:tmpNhi),5) = 2*hop.tpf(m)/1000*hop.results{m}{s}.hi(:,2);
            %remaining columns: means and stDevs in x/y and mean med_itrace
            tmpHi = zeros(tmpNhi,5);
            for n = 1:tmpNhi
                tmp_frames = hop.results{m}{s}.hi(n,1)+(0:hop.results{m}{s}.hi(n,2))-start_offset;
                tmp_XY = INPUT.XY{counter}(:,tmp_frames);
                tmpHi(n,1) = mean(tmp_XY(1,:));
                tmpHi(n,2) = std(tmp_XY(1,:));
                tmpHi(n,3) = mean(tmp_XY(2,:));
                tmpHi(n,4) = std(tmp_XY(2,:));
                tmpHi(n,5) = mean(INPUT.medI{counter}(tmp_frames));
            end
            allStats.hi(counterHi+(1:tmpNhi),6:10) = tmpHi;  
        end
        tmpNlo = size(hop.results{m}{s}.lo,1);
        if tmpNlo>0 
            allStats.lo(counterLo+(1:tmpNlo),1) = m;
            allStats.lo(counterLo+(1:tmpNlo),2) = hop.results{m}{s}.spotnum;
            allStats.lo(counterLo+(1:tmpNlo),3:4) = hop.results{m}{s}.lo;
            allStats.lo(counterLo+(1:tmpNlo),5) = 2*hop.tpf(m)/1000*hop.results{m}{s}.lo(:,2);
            tmpLo = zeros(tmpNlo,5);
            for n = 1:tmpNlo
                tmp_frames = hop.results{m}{s}.lo(n,1)+(0:hop.results{m}{s}.lo(n,2))-start_offset;
                tmp_XY = INPUT.XY{counter}(:,tmp_frames);
                tmpLo(n,1) = mean(tmp_XY(1,:));
                tmpLo(n,2) = std(tmp_XY(1,:));
                tmpLo(n,3) = mean(tmp_XY(2,:));
                tmpLo(n,4) = std(tmp_XY(2,:));
                tmpLo(n,5) = mean(INPUT.medI{counter}(tmp_frames));
            end
            allStats.lo(counterLo+(1:tmpNlo),6:10) = tmpLo;
        end           
        %update counters
        counterHi = counterHi + tmpNhi;
        counterLo = counterLo + tmpNlo;
        counter = counter + 1;
    end
end

% create output struct

output = struct('hop', hop, 'scatterStats', scatterStats, 'allStats', allStats);
output.XY = INPUT.XY;
output.medI = INPUT.medI;
            
    function [hi, lo] = get_hilo(traj, offset, ex_int)
        steps = find(diff(traj)~=0) + 1;
        if ~isempty(steps)
            % Start frames and lengths of states
            S = reshape(steps(1:end-1),length(steps)-1,1);
            L = reshape(steps(2:end)-steps(1:end-1),length(steps)-1,1);
            % Assign type of states
            updn = zeros(size(S));
            for i = 1:length(updn)
                updn(i) = sign(traj(steps(i))-traj(steps(i)-1));
            end
            % divide in hi and lo states
            hi = [S(updn==1)+offset L(updn==1)];
            lo = [S(updn==-1)+offset L(updn==-1)];
            if ~isempty(ex_int)
                hi = remove_excluded(hi, ex_int);
                lo = remove_excluded(lo, ex_int);
            end
        else
            hi = [];
            lo = [];
        end
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

    function RO = relative_occupancy(traj, K, ex_int, offset)
        if ~isempty(ex_int)
            tmp_ex = [];
            for i = 1:size(ex_int,1)
                tmp_ex = [tmp_ex (ex_int(i,1):ex_int(i,2))-offset];
            end
            tmp_ex(tmp_ex<0 | tmp_ex>length(traj))=[];
            if ~isempty(tmp_ex)
                traj(tmp_ex) = [];
            end
        end
        RO = sum(traj==K)/numel(traj);
    end
end





















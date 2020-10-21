function [ avg_frame ] = average_image(obj, varargin )
    % generate average image over a maximum of N_avg frames
    % start frame depends on input under 'mode' (default is frame #1)
    
    % parse input
    p = inputParser;
    
    addRequired(p, 'obj', @isobject);
    addOptional(p, 'N_avg', 100, @isnumeric);
    addOptional(p, 'start_frame', @isnumeric);
    addParameter(p, 'mode', 'first', @ischar);
    
    parse(p, obj, varargin{:});
    
    N_avg = p.Results.N_avg;
    if N_avg <= 0   % adjust to full movie length
        N_avg = length(obj.frames); 
    end
    
    switch p.Results.mode
        case 'first'
            N_max = N_avg;
            obj.initRead;
        case 'last'
            N_max = N_avg;
            obj.counter = length(obj.frames) - N_avg + 1;
        case 'from'
            N_max = min([N_avg length(obj.frames) - p.Results.start_frame + 1]);
            obj.counter = p.Results.start_frame;
    end    

    avg_frame = zeros(obj.sizeX, obj.sizeY);
    go_on = 1;
    N = 0;
    while go_on
        [movie, frames, go_on]  = obj.readNext;
        if N+length(frames) <= N_max
            avg_frame = avg_frame + sum(movie,3);
            N = N + length(frames);
        else
            avg_frame = avg_frame + sum(movie(:,:,1:(N_max-N)),3);
            go_on = 0;
            N = N_max;
        end
    end
    avg_frame = avg_frame ./ N; 

end
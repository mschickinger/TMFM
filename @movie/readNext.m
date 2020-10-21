function [tmp, frames_out, go_on] = readNext(obj)
    %reads the next N_read frames
    read_stop = obj.counter+obj.N_read-1; % last frame to read
    go_on = 1;
        if read_stop >= length(obj.frames)
            read_stop = length(obj.frames);
            go_on = 0;
        end

    tmp = zeros(obj.sizeX, obj.sizeY, read_stop-obj.counter+1);
    if isempty(obj.frames)
        frames_out = 1;
    else
        frames_out = obj.frames(obj.counter:read_stop);
    end 
    
    display(['Reading frame ' num2str(frames_out(1)) ' to ' num2str(frames_out(end))])
    
    for i=1:length(frames_out)
        tmp(:,:,i) = obj.readFrame(frames_out(i));
    end

    obj.counter = obj.counter + length(frames_out); 
end
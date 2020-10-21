function [itraces] = traces_movie_position(obj, positions, r_integrate)
% obtain intensity time traces by integrating specific ROIs in movie
    itraces = cell(0,1);
    go_on = 1;
    obj.initRead;
    while go_on
        [movie, frames, go_on]  = obj.readNext;
        itraces = append_traces_to_position(movie, itraces, frames, positions, r_integrate);
    end
end


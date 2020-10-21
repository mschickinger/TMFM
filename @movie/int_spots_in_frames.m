function [itraces] = int_spots_in_frames(obj, frame_nums, spots_pos, d_int)
% integrate intensities over specific ROIs in many frames (returns
        % spot intensity traces)
     itraces = cell(size(spots_pos,1),1);
     tmp = zeros(length(frame_nums),size(spots_pos,1));
     for j = 1:length(frame_nums)
         tmp(j,:) = obj.int_spots_in_frame(frame_nums(j), spots_pos, d_int);
     end
     for s = 1:size(itraces,1)
         itraces{s} = [reshape(frame_nums,length(frame_nums),1) zeros(length(frame_nums),3)];
         itraces{s}(:,2:3) = ones(length(frame_nums),2)*diag(spots_pos(s,:)) + obj.drift;
         itraces{s}(:,4) = tmp(:,s);
     end
end


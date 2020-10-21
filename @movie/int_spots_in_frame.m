function [ints] = int_spots_in_frame(obj, frame_idx, spots_pos, d_int)
 % integrate intensities over specific ROIs in one frame (suitable for parfor)
     ints = zeros(1,size(spots_pos,1));
     img = obj.readFrame(obj.frames(frame_idx));
     spots_pos = round(spots_pos + ones(size(spots_pos))*diag(obj.drift(frame_idx,:)));
     for i=1:size(spots_pos,1)
         sub_img = img(max(1,spots_pos(i,2)-d_int):min(obj.sizeY,spots_pos(i,2)+d_int),...
             max(1,spots_pos(i,1)-d_int):min(obj.sizeX,spots_pos(i,1)+d_int));
         ints(i) = sum(sub_img(:));
     end
end


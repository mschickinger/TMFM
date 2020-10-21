 function frames = getFrames(obj, sequence, first, last)
    %generate a list of images to be read from the movie
    frames = [];
    tmp = first:last;
    for i=1:length(tmp)
       if(  sequence(  mod(i-1,size(sequence,2))+1 )  )
          frames = [frames tmp(i) ];       
       end    
    end
 end

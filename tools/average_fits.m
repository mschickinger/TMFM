function avg_img = average_fits( filename, startframe, endframe )
%   Create average image of a substack of a .fits movie and output as
%   .tif and .png
info = fitsinfo(filename);
images = fitsread(filename, 'PixelRegion', {[1 info.PrimaryData.Size(1)], ...
    [1 info.PrimaryData.Size(2)],[startframe endframe]});
avg_img = mean(images,3);
end


function [h_min, p_out] = get_h_min(obj, r_find, N_img, varargin)
    %% determines peak-finding thresholds
    %   INPUTS: movie object
    %           'r_find' defines the spacing of the peaks
    %           'N_img' (can be omitted) gives number of images for average image
    %           'autosigma' (optional parameter) gives sigma multiplier
    %           value for automatic threshold
    %
    %   OUTPUTS:intensity threshold h_min
    %           positions of peaks with value higher than threshold
    % Example = get_h_min(obj, 4, 100, 'autosigma', 5)
    
    %% Parse input
    q = inputParser;
    
    addRequired(q, 'obj', @isobject);
    addRequired(q, 'r_find', @isnumeric);
    addOptional(q, 'N_img', @isnumeric);
    addParameter(q, 'autosigma', 0, @isnumeric);
    
    parse(q, obj, r_find, N_img, varargin{:});
    
    autosigma = q.Results.autosigma;
    
    %% Choose image for peak finding
    if exist('N_img', 'var') % use average frame
        img = obj.average_image(N_img); % average first N_img images
    else % use first frame if N_img is not specified
        if obj.input == 1 % tiff
            img = double(imread([obj.pname filesep obj.fnames{obj.frames(1)}]));
        else % fits
            img = fitsread([obj.pname filesep obj.fname{1}],  'Info', obj.info{1}, 'PixelRegion',{[1 obj.sizeX], [1 obj.sizeY], [obj.frames(1) obj.frames(1)] }); % read first frame                
        end
    end
    
    %% Find the peaks and choose threshold h_min
    p = find_peaks2d(img, r_find, 0, 0); % finding all possible peaks p has x, y, height, height-bg, I, I-I_bg
    
    img_std = std(img(:));
   
    if autosigma % true if autosigma > 0
        h_min = autosigma*img_std;
    else
        p_7std = p(p(:,4)>=7*img_std, :); % this is just an estimate, #of peaks found may vary since peak_find algorithm return height as int not double
        p_5std = p(p(:,4)>=5*img_std, :);
        p_3std = p(p(:,4)>=3*img_std, :);

        % plot
        close all
        figure('units','normalized','outerposition',[0 0 1 1])
        subplot(1, 2, 1)
        imagesc(img), colorbar, axis image, colormap gray, hold on
        if size(p_3std,1)>0
            h(1) = plot(p_3std(:,1)+1, p_3std(:,2)+1, 'ro');
            if size(p_5std,1)>0
                h(2) = plot(p_5std(:,1)+1, p_5std(:,2)+1, 'go');
            end
            if size(p_7std,1)>0
                h(3) = plot(p_7std(:,1)+1, p_7std(:,2)+1, 'bo');
            end
            legend(h, {['3\sigma = ' num2str(round(3*img_std)) ], ['5\sigma = ' num2str(round(5*img_std)) ], ['7\sigma = ' num2str(round(7*img_std)) ]})
        else
            h = plot(p(p(:,4)>=.5*img_std,1)+1, p(p(:,4)>=.5*img_std,2)+1,'ro');
            legend(h, {['\frac{\sigma}{2} = ' num2str(round(.5*img_std))]});
        end   

        subplot(1, 2, 2)
        xhist = min(p(:,4)):5:max(p(:,4));
        n = hist(p(:,4), xhist);
        semilogy(xhist, sum(n)-cumsum(n)), hold on
        h(1) = vline(3*img_std, 'r');
        h(2) = vline(5*img_std, 'g');
        h(3) = vline(7*img_std, 'b');
        legend(h, {['3\sigma = ' num2str(round(3*img_std)) ], ['5\sigma = ' num2str(round(5*img_std)) ], ['7\sigma = ' num2str(round(7*img_std)) ]})
        set(gca, 'XLim', [0 xhist(end)])
        xlabel('Minimal height'), ylabel('# of peaks found')
        axis square

        options.WindowStyle='normal';
        prompt={'Enter min height (default=5*sigma):'};
        def={num2str(round(5*img_std))};
        threshold =inputdlg(prompt, strcat('Enter threshold:'), 1,   def, options);  % leave out
        h_min =str2double(threshold(1)); % leave out
        pause(0.1) % leave out
        close all
    end
    
    obj.h_min = h_min;
    p_out = p(p(:,4)>=h_min,:);
end
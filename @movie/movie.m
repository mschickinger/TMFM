classdef movie < handle
    %TRACE_SELECTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        frames; % all images to be read
        N_read = 300; % number of images to read in one portion
        counter = 1; % internal counter for reading the movie
        
        sequence; % sequence to be read, e.g. 101010
        first; % first image to be read
        last; % last image to be read
        
        pname; % pathname of file location
        fname; % filename of file
        
        sizeX; % number of pixel in X-dimension
        sizeY; % number of pixel in Y-dimension
        mov_length; % number of frames in the whole movie
        
        info; % fits info
        h_min; % minimal height for peak finding
        
        input; % 0=fits, 1=tiff-stack
        fnames; % cell with all filenames, only for tiff-stack

        N_frame_per_fits; % stores the number of frames in one fits file
        add_upon_reading; % value added to frame after fitsread. Depends on .fits intercept
        
        drift; % stores displacement in x and y over whole movie through drift.
    end
    
    methods
        %constructor
        function obj = movie(pname, fname, first, last, sequence) % fname is the filename of the first fits-file
            obj.sequence = sequence;
            obj.pname = pname;
            obj.fname = cell(1,1);
            obj.fname{1} = fname;
            
            if strcmp(fname(end-2:end), 'tif')
                obj.input = 1; % read tiff data
            else
                obj.input = 0; % read fits data
            end
            
            if obj.input == 1 % tiff-stack
                tmp = dir([pname filesep '*.tif']);
                obj.fnames = {tmp.name};
                obj.info = 'Tiff-Stack info';
                obj.sizeX = size(imread([pname filesep obj.fnames{1}]),2);
                obj.sizeY = size(imread([pname filesep obj.fnames{1}]),1);
                obj.mov_length = length(obj.fnames);
            else % fits
                obj.info = cell(1,1);
                obj.info{1} = fitsinfo([obj.pname filesep obj.fname{1}]);
                
                obj.N_frame_per_fits = 4095; %obj.info{1}.PrimaryData.Size(3);
                
                if obj.info{1}.PrimaryData.Intercept==0
                    tmp = fitsread([pname filesep fname],'PixelRegion',{[1 obj.sizeX],[1 obj.sizeY],[1 1]});
                    if max(tmp(:)) < -2^8 && min(tmp(:))>-2^16;
                        obj.add_upon_reading = 2^16;
                    else
                        obj.add_upon_reading = 0;
                    end
                else
                    obj.add_upon_reading = 0;
                end
                
                obj.sizeX = obj.info{1}.PrimaryData.Size(1); 
                obj.sizeY = obj.info{1}.PrimaryData.Size(2);
                
                tmp = dir([pname filesep fname(1:end-4) '_X*.fits']); % change * to wildcard for 1-2 character/integers
                [~,idx] = sort([tmp.datenum]);
                tmp = tmp(idx);
                for i=1:length(tmp)
                    obj.fname = [obj.fname tmp(i).name];
                    obj.info = [obj.info fitsinfo([obj.pname filesep tmp(i).name])];
                end
                
                f_tot = 0; % calculate total number of frames
                for i=1:length(obj.fname)
                    switch length(obj.info{i}.PrimaryData.Size)
                        case 3
                            f_tot = f_tot + obj.info{i}.PrimaryData.Size(3); 
                        case 2
                            f_tot= f_tot + 1;
                    end           
                end
                obj.mov_length = f_tot;
            end
            obj.first = first;
            if last == -1
                obj.last = obj.mov_length;
            else
                obj.last = last;
            end
            obj.frames = obj.getFrames(sequence, obj.first, obj.last);
            disp('Movie class created.')
        end
        
        % generate a list of images to be read from the movie
        frames = getFrames(obj, sequence, first, last)
        
        % reads one frame
        [img] = readFrame(obj,framenumber)
        
        % reads the next N_read frames
        [tmp, frames_out, go_on] = readNext(obj)
        
        % initialize the counter
        initRead(obj)
        
        % generate average image, starting from first frame until N_max
        [avg_frame] = average_image(obj, varargin)
        
        % determine peak-finding thresholds

        [h_min, p_out] = get_h_min(obj, r_find, N_img, varargin)
        
        % obtain intensity time traces by integrating specific ROIs in movie
        itraces = traces_movie_position(obj, positions, r_integrate)
        
        % integrate intensities over specific ROIs in one frame (suitable for parfor)
        ints = int_spots_in_frame(obj, frame_idx, spots_pos, d_int)
        
        % integrate intensities over specific ROIs in many frames (returns
        % spot intensity traces)
        itraces = int_spots_in_frames(obj, frame_nums, spots_pos, d_int)
    
    end
      
end        
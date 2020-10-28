function [  ] = stats_to_igor_nost( stats, ID, path )
% Writes a txt file that can be imported to igor pro for a scatter plot
    
if nargin == 2
    path = cd;
end
wave_names = {[ID '_mTu_nost'],[ID '_mTb_nost'],[ID '_SEMu_nost'],[ID '_SEMb_nost']};
fullpath = [path filesep ID '_Scatter_nonstick.txt'];
    if size(stats,2)== length(wave_names)
        % write header/wave names
        fileID=fopen(fullpath, 'w'); %open file to write
        for i=1:size(stats,2)            %write wavenames at each column header
            fprintf(fileID, [wave_names{i} '\t']);
        end
        fprintf(fileID,'\n');
        fclose(fileID);

        % append data
        dlmwrite(fullpath, stats, 'delimiter', '\t','-append')
    else
        disp('Error: columns not equal to number of wave names.')
    end

end



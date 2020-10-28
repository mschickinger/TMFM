function [allD, allDmax] = get_allDmax(densities, discard)%, maxlBelow)
    allD = cell(2,1);
    allDmax = cell(2,1);
    for k = 1:2
        for isp = setdiff(1:size(densities,1),discard)
            if ~isempty(densities{isp,k})
                allD{k} = [allD{k};densities{isp,k}(~isnan(densities{isp,k}(:,1)),1)];
                allDmax{k} = [allDmax{k};densities{isp,k}(~isnan(densities{isp,k}(:,2)),2)];
            end
        end
    end

end


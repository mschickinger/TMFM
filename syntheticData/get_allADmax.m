function [allD, allDmax, allA, allAmax, allLmax] = get_allADmax(densities, discard, areas)%, maxlBelow)
    allD = cell(2,1);
    allDmax = cell(2,1);
    allA = cell(2,1);
    allAmax = cell(2,1);
    %allLmax = cell(2,1);
    for k = 1:2
        for isp = setdiff(1:size(densities,1),discard)
            if ~isempty(densities{isp,k})
                allD{k} = [allD{k};densities{isp,k}(~isnan(densities{isp,k}(:,1)),1)];
                allDmax{k} = [allDmax{k};densities{isp,k}(~isnan(densities{isp,k}(:,2)),2)];
                allA{k} = [allA{k};areas{isp,k}(~isnan(areas{isp,k}(:,1)),1)];
                allAmax{k} = [allAmax{k};areas{isp,k}(~isnan(areas{isp,k}(:,2)),2)];
                %allLmax{k} = [allLmax{k};maxlBelow{isp,k}(~isnan(maxlBelow{isp,k}))];
            end
        end
    end

end


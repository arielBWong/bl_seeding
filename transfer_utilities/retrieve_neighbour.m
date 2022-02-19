function [neighbour_optxu, neighbour_optxl, close_id] = retrieve_neighbour(xu, lower_searchdata_cly, archive_xu, archive_xl)
% this function find a neighbour xu from all neigbour archive
% this neighbour has a lower level search record of global search (EA)
% input:
% xu is current variable
% lower_searchdata_cly is a cell array recording all data evaluated on the
% lower leve for every xu in archive_xu
% archive_xu, bilevel search archive xu
% archive_xl, bilevel search archive xl
% output:
% neigbhour and its matching xl
% ------------------

global lower_decisionSwitch
dist = pdist2(xu, archive_xu); % this xu is upper level new infill xu, not added into archive_xu
[~, idx] = sort(dist);
for i = 1: length(idx)
    % if size(lower_searchdata_cly{idx(i)}, 1) > 50
    if lower_decisionSwitch(idx(i)) == 0      
        close_id = idx(i);
        break
    end
end

if size(lower_searchdata_cly{close_id}, 1) < 50
    a = 0;
end

neighbour_optxu = archive_xu(close_id, :);
neighbour_optxl = archive_xl(close_id, :);

end
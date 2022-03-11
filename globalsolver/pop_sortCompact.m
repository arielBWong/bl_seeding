function [pop] = pop_sortCompact(pop)
% auxiliary function
% this function sort evoluation popluation w.r.t.
% number of objectives
% constraints
% nd sort is not compatible with single objective problems
%-----------------------------------------------------

numcon = size(pop.C, 2);
numobj = size(pop.F, 2);

if numcon == 0
    if numobj > 1                                                  % mo problem
        [fronts, ids, ~] = nd_sort(pop.F, (1:size(pop.F, 1))');
    else                                                                        % so problem
        [~, ids] = sort(pop.F);                                          % acending sort/minimization
    end
    
    %------------
    
    pop.F = pop.F(ids, :);
    pop.X = pop.X(ids, :);
    
    if ~isempty(pop.A)
        pop.A =  pop.A(ids,:);
    end
    
    if ~isempty(pop.Mdl)
        pop.Mdl = pop.Mdl(ids);   % cell array
    end
    
    if ~isempty(pop.trgdata)
        pop.trgdata = pop.trgdata(ids);   % cell array
    end
    
    if ~isempty(pop.switch_lls)
        pop.switch_lls = pop.switch_lls(ids);   % cell array
    end
    
    
    
    
end

if numcon ~= 0
    error('not suitable for constraint problems');
end

end

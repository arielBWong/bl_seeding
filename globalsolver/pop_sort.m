function [sf, sx, sc, sa] = pop_sort(f, x, c, a)
% auxiliary function
% this function sort evoluation popluation w.r.t.
% number of objectives
% constraints
% nd sort is not compatible with single objective problems
% input
%               f                                                       %  objective population
%               x                                                       %  design variable population
%               c                                                       %  constraints population
%               a
%              
%               
% output
%               sf                                                      % sorted objective population
%               sx                                                      % sorted design variable population
%               sc                                                      % sorted constraints population
%               sa
%-----------------------------------------------------

numcon = size(c, 2);
numobj = size(f, 2);

if numcon == 0                                                  % unconstraint problem
    sc = [];
    if numobj > 1                                               % mo problem
        [~, ids, ~] = nd_sort(f, (1:size(f, 1))');
    else                                                        % so problem
        [~, ids] = sort(f);                                     % acending sort/minimization
    end
    sf = f(ids,:);
    sx = x(ids, :);
    % process addon
    if ~isempty(a)
        sa =  a(ids,:);
    else
        sa = [];
    end
    return;
end

sa = []; % for cases sa is not used
if numcon>0
    % Ordering considers constraints
    c(c<=0) = 0;
    fy_ind = sum(c, 2) == 0;                 % feasibility index
    cv_ind =~fy_ind;                                            % constraint violation index
    
    
    % seperate feasible and infeasible
    % sort two subset seperately
    fy_F = f(fy_ind, :); cv_F = f(cv_ind, :);                % operation should be valid when no feasible solutions
    fy_C = c(fy_ind, :); cv_C = c(cv_ind, :);
    fy_X = x(fy_ind, :); cv_X = x(cv_ind, :);
    
    % process addon
    if ~isempty(a)
        fy_A = a(fy_ind, :); cv_A = a(cv_ind, :); 
    end
    
    
    % sort feasible
    if numobj>1
        [~, ids, ~] = nd_sort(fy_F, (1: size(fy_F, 1))');   % reason to do this is, nd_sort.m is not compatible with so
    else
        [~, ids] = sort(fy_F);
    end
    fy_F = fy_F(ids, :); fy_C = fy_C(ids, :); fy_X = fy_X(ids, :);
    if ~isempty(a)
        fy_A = fy_A(ids, :);
    end
    
    % sort infeasible
    if numcon > 1
        sum_cv = sum(cv_C, 2);
    else
        sum_cv = cv_C;
    end
    [~, idc] = sort(sum_cv); 
    cv_F = cv_F(idc, :); cv_C = cv_C(idc, :); cv_X = cv_X(idc, :);
    if ~isempty(a)
        cv_A = cv_A(idc, :);
    end
    
    % replace unsorted each fields of pop
    sf = [fy_F; cv_F]; sc= [fy_C; cv_C]; sx = [fy_X; cv_X]; 
    if ~isemtpy(a)
        sa = [fy_A; cv_A];
    end
end
end

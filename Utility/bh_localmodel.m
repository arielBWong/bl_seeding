function [localx, localf, localc, krg, krgc, arc_obj, arc_c] = bh_localmodel(x, f, c, varargin)
% create a good local region around feasible min f
% return local regions' data, model built from training data
%----------------------------

prob            = varargin{1};

% feasibility identification
if ~isempty(c)
    % change f to order
    [sx, sf, sc, s] ...
               = sort_archive(x, f, c);
else
    sx = x; sf = f; sc = c; s = f;
end



% start to locate a local region
prim            =  bump_detection(sx, s);
nboxes          =  length(prim.boxes);
nvar            =  size(sx, 2);
surrogate_minsize ...
                = 2 * (nvar + 1);


upper_stack     = [];
lower_stack     = [];
localsurrogate  = false;

% identify minf
[~, id_min]     = min(s, [], 1);
minx            = sx(id_min, :);

for ii = 1: nboxes -1
    [bump_lb, bump_ub] = boxboundary(prim, ii, prob);
    upper_stack        = [upper_stack; bump_ub];
    lower_stack        = [lower_stack; bump_lb];
    
    % adjust upper bound and lower bound
    upper_bound        = max(upper_stack, [], 1);
    lower_bound        = min(lower_stack, [], 1);
    
    % count enclosed points
    [localx, localf, localc] ...
                       = pickup_localxy(upper_bound, lower_bound, sx, sf, sc);
                   
    % identify region quality
    [localx, localf, localc, train_idx] ...
                       = data_prepare(localx, localf, localc);
    
  
    
    if size(localx, 1) >= surrogate_minsize                  % minimum number of local search point
        if inboxcheck(minx, upper_bound, lower_bound)       % fmin included
            localsurrogate     = true;
            % kbs = size(localx, 1);
            % fprintf('local search size %d\n', kbs);
            break;
        end
    end
end



if localsurrogate
    [krg, krgc,arc_obj, arc_c] = surrogate_build(localx, localf, localc, train_idx);  
else
    krg = {}; % signal for outside program to escape local search
    localx = []; localf=[]; localc=[]; krgc={}; arc_obj=[];arc_c=[];
end
end


function [inbox_flag] = inboxcheck(x,  bu, lb)
    inbox_flag = false;
    nx = size(x,2); 
   if sum(bu - x >= 0 ) == nx && sum(x - lb >= 0)
       inbox_flag = true;
   end
    
end


function [sx, sf, sc, s] = sort_archive(x, f, c)
% sort archive (x, f, c ) add one more index s for replacing objecive
numobj  = size(f, 2);
numcon  = size(c, 2);
% Ordering considers constraints
c(c<=0) = 0;
fy_ind  = sum(c, 2) == 0;                                   % feasibility index
cv_ind  = ~fy_ind;                                          % constraint violation index

    
% seperate feasible and infeasible
% sort two subset seperately
fy_F = f(fy_ind, :); cv_F = f(cv_ind, :);                   % operation should be valid when no feasible solutions
fy_C = c(fy_ind, :); cv_C = c(cv_ind, :);
fy_X = x(fy_ind, :); cv_X = x(cv_ind, :);

% sort feasible
if numobj>1
    [~, ids, ~] = nd_sort(fy_F, (1: size(fy_F, 1))');      % reason to do this is, nd_sort.m is not compatible with so
else
    [~, ids]    = sort(fy_F);
end
fy_F = fy_F(ids, :); fy_C = fy_C(ids, :); fy_X = fy_X(ids, :);
    
% sort infeasible
if numcon > 1
    sum_cv     = sum(cv_C, 2);
else
    sum_cv     = cv_C;
end
[~, idc]       = sort(sum_cv);
cv_F = cv_F(idc, :); cv_C = cv_C(idc, :); cv_X = cv_X(idc, :);

% replace unsorted each fields of pop
sf = [fy_F; cv_F]; sc = [fy_C; cv_C]; sx = [fy_X; cv_X];
s = 1: size(sf, 1);
s = s';

end


function [localx, localy, localc] = pickup_localxy(upper_bound, lower_bound, xx, yy, cc)
localx = [];
localy = [];
localc = [];

numx   = size(xx, 1);
nx     = size(xx, 2);
for ii = 1: numx
    x = xx(ii, :);
    y = yy(ii, :);
    if ~isempty(cc)
        c = cc(ii, :); end
    
    % with in boundary
    if sum(upper_bound - x >= 0 ) == nx && sum(x - lower_bound >= 0) == nx
        localx = [localx; x];
        localy = [localy; y];
        if  ~isempty(cc)
            localc = [localc; c];
        end
    end
end
end


function[bump_lb, bump_ub] = boxboundary(prim, ii, prob)
mentioned_var = prim.boxes{ii}.vars;
bump_lb = prob.xl_bl;
bump_ub = prob.xl_bu;

nb = length(mentioned_var);
if ~isempty(mentioned_var)
    for jj = 1:nb % varible indicator
        if ~isnan(prim.boxes{ii}.min(jj))
            bump_lb(mentioned_var(jj)) =  prim.boxes{ii}.min(jj);
        end
        
        if ~isnan(prim.boxes{ii}.max(jj))
            bump_ub(mentioned_var(jj)) =  prim.boxes{ii}.max(jj);
        end 
    end
end
end



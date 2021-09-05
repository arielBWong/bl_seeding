
function [best_x] = localsurrgoate_search(krg, krgc, local_ub, local_lb, arc_obj, arc_c, num_pop, num_gen, varargin)
% search on surrogate 
%

funh_obj = @(x)surrogate_predict(x, krg, arc_obj);
funh_con = @(x)surrogate_predict(x, krgc, arc_c);

param.gen        = num_gen;
param.popsize    = num_pop;


l_nvar           = size(local_ub, 2);

lb               = varargin{1};
ub               = varargin{2};

expand           = 0.05 * abs(local_ub - local_lb);
local_ub         = min(local_ub + expand, ub);
local_lb         = max(local_lb - expand, lb);

visual           = [];
% expand local search to 
[~,~,~, archive] = gsolver(funh_obj, l_nvar,  local_lb, local_ub, [], funh_con, param);
best_x           = archive.pop_last.X(1, :);
end
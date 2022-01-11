function [starting_xl]= cokriging_decisionmaking(co_mdl, prob, cokrg_lb, cokrg_ub, varargin)
% global search on cokrging to decide a starting point
p = inputParser;
addRequired(p, 'co_mdl');
addRequired(p, 'prob');
addRequired(p, 'cokrg_lb');
addRequired(p, 'cokrg_ub');
addParameter(p, 'lx', []);
addParameter(p, 'lf', []);                           % probing archive
addParameter(p, 'hx', []);                         % conduct seeding controller
addParameter(p, 'hf', []);                                 % additional restart times
addParameter(p, 'boundary_updatesize', 0);
addParameter(p, 'xu', []);
addParameter(p, 'extend', false);


parse(p, co_mdl, prob, cokrg_lb, cokrg_ub, varargin{:});
% -----------------------------------

co_mdl = p.Results.co_mdl;
prob = p.Results.prob;
cokrg_lb = p.Results.cokrg_lb;
cokrg_ub = p.Results.cokrg_ub;
cheap_x = p.Results.lx;
cheap_f = p.Results.lf;
expensive_x = p.Results.hx;
expensive_f = p.Results.hf;
boundary_updatesize = p.Results.boundary_updatesize;
xu = p.Results.xu;
extend = p.Results.extend;
%-----------------------------------------


funh_obj = @(x)current_landscapeF(x, co_mdl);
funh_con = @(x)noconstraint(x);

num_xvar = prob.n_lvar;

ea_param.gen = 100; 
ea_param.popsize = 100;

[starting_xl, ~, ~, ~, ~] = gsolver(funh_obj, num_xvar, cokrg_lb, cokrg_ub, [], funh_con, ea_param);


[onbound_check, new_lb, new_ub] = onbound(starting_xl, cokrg_lb, cokrg_ub, prob);
if onbound_check && extend
    new_xl        = lhsdesign(boundary_updatesize, prob.n_lvar, 'criterion','maximin','iterations',1000);
    new_xl        = repmat(new_lb, half_global, 1) ...
                        + repmat((new_ub - new_lb), boundary_updatesize, 1) .* new_xl;
    xu_ext           = repmat(xu, boundary_updatesize, 1);
    new_fl        = prob.evaluate_l(xu_ext, new_xl);
    
    expensive_x = [expensive_x; new_xl];
    expensive_f = [expensive_f; new_fl];
    
    [~, ia, ~]  = unique(expensive_x, 'rows');    % check unique
    expensive_x = expensive_x(ia, :);
    expensive_f = expensive_f(ia, :);
    [expensive_x, expensive_f] = keepdistance_rawdata(expensive_x, expensive_f, prob.xl_bl, prob.xl_bu);
    
    
    x_trg{1,:} = cheap_x;                    % cheap
    y_trg{1,:} = cheap_f;
    x_trg{2,:} = expensive_x;
    y_trg{2,:} = expensive_f;
    
    
    co_mdl = oodacefit(x_trg, y_trg);
    funh_obj = @(x)current_landscapeF(x, co_mdl);
    funh_con = @(x)noconstraint(x);
    [starting_xl, ~, ~, ~, ~] = gsolver(funh_obj, num_xvar, new_lb, new_ub, [], funh_con, ea_param);
    
end

end

function [f] = current_landscapeF(x, mdl)

[f,~] = mdl.predict(x);
end

function c = noconstraint(x)
c = [];
end

function [flag, new_lb, new_ub] = onbound(x, lb, ub, prob)
 

flag = false;
check_ub = ub - x;
id_ub = check_ub < 1d-6;
if any(id_ub > 0)
    flag = true;
end

check_lb = x - lb;
id_lb = check_lb < 1e-6;
if any(id_lb > 0)
    flag = true;
end

new_lb = prob.xl_bl .*  id_lb + lb;
new_ub = prob.xl_bu .*  id_ub + ub;

end



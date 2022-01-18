function [starting_xl, onbound_check]= cokriging_decisionmaking(co_mdl, prob, cokrg_lb, cokrg_ub, varargin)
% global search on cokrging to decide a starting point
p = inputParser;
addRequired(p, 'co_mdl');
addRequired(p, 'prob');
addRequired(p, 'cokrg_lb');
addRequired(p, 'cokrg_ub');

parse(p, co_mdl, prob, cokrg_lb, cokrg_ub, varargin{:});
% -----------------------------------
co_mdl = p.Results.co_mdl;
prob = p.Results.prob;
cokrg_lb = p.Results.cokrg_lb;
cokrg_ub = p.Results.cokrg_ub;
%-----------------------------------------


funh_obj = @(x)current_landscapeF(x, co_mdl);
funh_con = @(x)noconstraint(x);

num_xvar = prob.n_lvar;

ea_param.gen = 100; 
ea_param.popsize = 100;

[starting_xl, ~, ~, ~, ~] = gsolver(funh_obj, num_xvar, cokrg_lb, cokrg_ub, [], funh_con, ea_param);
[onbound_check] = onbound(starting_xl, cokrg_lb, cokrg_ub, prob);

% boundary process is: if hit boundary, re-evaluate
% re-create starting point 
% if still hit boundary, check problem boundary
% if on problem boudnary, then accept it



end

function [f] = current_landscapeF(x, mdl)

[f,~] = mdl.predict(x);
end

function c = noconstraint(x)
c = [];
end

function [flag] = onbound(x, lb, ub, prob)
% this function  check in normalized space
% whether x is on boundary

flag_up = false;
flag_down = true;

x_norm = (x - prob.xl_bl)./(prob.xl_bu - prob.xl_bl);
lb_norm = (lb - prob.xl_bl) ./(prob.xl_bu - prob.xl_bl);
ub_norm = (ub-prob.xl_bl) ./(prob.xl_bu - prob.xl_bl);


check_ub = ub_norm - x_norm;
id_ub = check_ub < 1e-6;
if any(id_ub > 0)
   flag_up = true; 
end

check_lb = x_norm - lb_norm;
id_lb = check_lb < 1e-6;
if any(id_lb > 0)
    flag_down = true;
end

flag = flag_up|flag_down;
end



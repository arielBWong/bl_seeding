function [starting_xl, onbound_check]= cokriging_decisionmaking(co_mdl, prob, cokrg_lb, cokrg_ub, varargin)
% global search on cokrging to decide a starting point
p = inputParser;
addRequired(p, 'co_mdl');
addRequired(p, 'prob');
addRequired(p, 'cokrg_lb');
addRequired(p, 'cokrg_ub');
addRequired(p, 'bd_funh');

parse(p, co_mdl, prob, cokrg_lb, cokrg_ub, varargin{:});
% -----------------------------------
co_mdl = p.Results.co_mdl;
prob = p.Results.prob;
cokrg_lb = p.Results.cokrg_lb;
cokrg_ub = p.Results.cokrg_ub;
bd_funh = p.Results.bd_funh;
%-----------------------------------------


% funh_obj = @(x)current_landscapeF(x, co_mdl);
funh_obj = @(x)current_landscapeF(x, co_mdl, prob);
funh_con = @(x)noconstraint(x);

num_xvar = prob.n_lvar;

ea_param.gen = 100; 
ea_param.popsize = 100;

[starting_xl, ~, ~, ~, ~] = gsolver(funh_obj, num_xvar, cokrg_lb, cokrg_ub, [], funh_con, ea_param);
[onbound_check] = bd_funh(starting_xl);


end





function [f] = current_landscapeF(x, mdl, prob)
% x = (x - prob.xl_bl) ./ (prob.xl_bu - prob.xl_bl);
[f, ~] = mdl.predict(x);
end

function c = noconstraint(x)
c = [];
end



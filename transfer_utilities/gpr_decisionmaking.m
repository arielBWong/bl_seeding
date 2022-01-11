function [starting_xl]= gpr_decisionmaking(mdl, prob, xu)


funh_obj = @(x)gpr_prediction(x, mdl, xu);
funh_con = @(x)noconstraint(x);

num_xvar = prob.n_lvar;

ea_param.gen = 100; 
ea_param.popsize = 100;

[starting_xl, ~, ~, ~, ~] = gsolver(funh_obj, num_xvar, prob.xl_bl, prob.xl_bu, [], funh_con, ea_param);

end

function [f] = gpr_prediction(x, mdl, xu)
n = size(x, 1);
xu = repmat(xu, n, 1);

fullx = [xu, x];
[f,~] = predict(mdl, fullx);
end

function c = noconstraint(x)
c = [];
end

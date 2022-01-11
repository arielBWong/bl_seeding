function [starting_xl,  mdl_delta]= correctionFunc_decisionmaking(prob_xl, fl_current, fl_close, mdl_close, prob)
% this function uses surrogate to create global search 
delta_d = fl_current - fl_close;

param.GPR_type = 2;
param.no_trials = 1;
mdl_delta = Train_GPR(prob_xl, delta_d, param);

funh_obj = @(x)current_landscapeF(x, mdl_delta, mdl_close);
funh_con = @(x)noconstraint(x);

num_xvar = prob.n_lvar;

ea_param.gen = 100; 
ea_param.popsize = 100;

[starting_xl, ~, ~, ~, ~] = gsolver(funh_obj, num_xvar, prob.xl_bl, prob.xl_bu, [], funh_con, ea_param);

end


function [f] = current_landscapeF(x, mdl_delta, mdl_closeland)

param.GPR_type   = 2;
param.no_trials  = 1;

[delta_f, delta_mu] = Predict_Objective(mdl_delta, x, param);
[close_f, close_mu] = Predict_Objective(mdl_closeland, x, param);

f = close_f + delta_f;
end


function c = noconstraint(x)
c = [];
end

% Predicting using a Kriging model
function [y_prd, s_prd] = Predict_Objective(gprMdl,x_tst,param)
y_prd=[]; s_prd=[];
for i=1:length(gprMdl)
    if(param.GPR_type==1)
        [tmp_mu,tmp_sigma] = predict(gprMdl{i},x_tst);
    else
        [tmp_mu,tmp_sigma] = Predict_DACE(gprMdl{i},x_tst,i);
    end
  
    y_prd = [y_prd tmp_mu];
    s_prd = [s_prd tmp_sigma];
end
end

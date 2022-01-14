function[match_xl, lower_evalcount, x_visited, f_visited] = localsearch_withTransferredSolution(prob, xu, tf_xl, maxFE, ...
    lower_evalcount, varargin)
% this function  takes  a starting point and conduct local search on true
% evaluation, it returns local search result xl, total function evaluation
% numbers and solutions evaluated (TO NOTE: this two numbers are different due to fmincon interal mechanism)
% 
%
%----------------------------
p = inputParser;
addRequired(p, 'prob');   % problem
addRequired(p, 'xu');   % xu
addRequired(p, 'tf_xl');  % transferred solution
addRequired(p, 'maxFE');  % FE budget
addRequired(p, 'lower_evalcount');  % Number of FE used
parse(p, prob, xu, tf_xl, maxFE, lower_evalcount, varargin{:});
% -----------------------------------
prob = p.Results.prob;
xu = p.Results.xu;
tf_xl = p.Results.tf_xl;
maxFE = p.Results.maxFE;
lower_evalcount = p.Results.lower_evalcount;
%-------------------------

funh_obj = @(x)objective_func(prob, xu, x);
funh_con = @(x)constraint_func(prob, xu, x);

[tf_fl, tf_fc] = prob.evaluate_l(xu, tf_xl);
[match_xl, ~, history, output] = lowerlevel_fmincon(tf_xl, maxFE, prob.xl_bl, prob.xl_bu,funh_obj, funh_con);

[match_fl, match_cl] = prob.evaluate_l(xu, match_xl);  % additional lazy step, can be extracted from local search results
additional_searchxl = [tf_xl; match_xl];   
additional_searchfl = [tf_fl; match_fl];
additional_searchcl = [tf_fc; match_cl];

[match_xl, ~, ~, ~, ~] =  localsolver_startselection(additional_searchxl, additional_searchfl, additional_searchcl);

% Collect solutions from all evaluated solutions
% local search might not have as many as counted
x_visited = history.x;
f_visited = history.fval;
lower_evalcount = lower_evalcount + output.funcCount;

end


function f = objective_func(prob, xu, xl)
f = prob.evaluate_l(xu, xl);
end

function [c, ce] = constraint_func(prob, xu, xl)
c = [];
ce = [];
end
function[match_xl, lower_eval, lowerlocal_record, x_visited, f_visited] = localsearch_withTransferredSolution(prob, xu, tf_xl, maxFE, ...
    lower_eval, lowerlocal_record, varargin)
% 
%
%

p = inputParser;
addRequired(p, 'prob'); % problem
addRequired(p, 'xu');   % xu
addRequired(p, 'tf_xl'); % transferred solution
addRequired(p, 'maxFE'); % FE budget
addParameter(p, 'lower_eval', 0);  % Number of FE used
addParameter(p, 'lowerlocal_record', []); % test purpose for checking whether local search is success or not

addParameter(p, 'close_optxl',  []);
addParameter(p, 'close_optxu',  []);
addParameter(p, 'close_optxlcokring',[]);
addParameter(p, 'cheap_x', []);
addParameter(p, 'cheap_f', []);
addParameter(p, 'expensive_x', []);
addParameter(p, 'expensive_f', []);
addParameter(p, 'co_mdl', []);
addParameter(p, 'close_archive', [])


parse(p, prob, xu, tf_xl, maxFE, lower_eval, lowerlocal_record, varargin{:});
% -----------------------------------

prob = p.Results.prob;
xu = p.Results.xu;
tf_xl = p.Results.tf_xl;
maxFE = p.Results.maxFE;
lower_eval = p.Results.lower_eval;
lowerlocal_record = p.Results.lowerlocal_record;

close_optxl = p.Results.close_optxl;
close_optxu= p.Results.close_optxu;
close_optxlcokring = p.Results.close_optxlcokring;
cheap_x = p.Results.cheap_x;
cheap_f = p.Results.cheap_f;
expensive_x = p.Results.expensive_x;
expensive_f = p.Results.expensive_f;
co_mdl = p.Results.co_mdl;
close_archive = p.Results.close_archive;
%-------------------------


funh_obj = @(x)objective_func(prob, xu, x);
funh_con = @(x)constraint_func(prob, xu, x);

[tf_fl, tf_fc] = prob.evaluate_l(xu, tf_xl);
[match_xl, ~, history, output] = lowerlevel_fmincon(tf_xl, maxFE, prob.xl_bl, prob.xl_bu,funh_obj, funh_con);

lower_eval = lower_eval + output.funcCount;


[match_fl, match_cl] = prob.evaluate_l(xu, match_xl);  % additional lazy step, can be extracted from local search results
additional_searchxl = [tf_xl; match_xl];               % variable name from copy paste
additional_searchfl = [tf_fl; match_fl];
additional_searchcl = [tf_fc; match_cl];

[match_xl, match_fl, ~, flag, ~] =  localsolver_startselection(additional_searchxl, additional_searchfl, additional_searchcl);


% for test purpose, whether local search success
ideal_xl = prob.get_xlprime(xu);
[ideal_fl, ~] = prob.evaluate_l(xu, ideal_xl);
if abs(ideal_fl - match_fl) > 1
    lowerlocal_record = [lowerlocal_record; 0]; % local search fail
    
    % plot2d_withCokring(xu, prob.xl_bl, prob.xl_bu, match_xl, close_optxl, prob, close_optxu, co_mdl,  close_optxlcokring, ...
                       %cheap_x, cheap_f, expensive_x, expensive_f, close_archive);
else
    lowerlocal_record = [lowerlocal_record; 1]; % local search success
    
     %plot2d_withCokring(xu, prob.xl_bl, prob.xl_bu, match_xl, close_optxl, prob, close_optxu, co_mdl,  close_optxlcokring, ...
     %                    cheap_x, cheap_f, expensive_x, expensive_f, close_archive);
end



% Collect solutions from all evaluated solutions
% local search might not have as many as counted
x_visited = history.x;
f_visited = history.fval;


% if correctionFunc_decision   %  else mdl='position_holder' set at the beginning
%     mdl = lower_surrogateCreate(aftersearch_collectx, aftersearch_collectf, prob.xl_bl, prob.xl_bu);
% end


end


function f = objective_func(prob, xu, xl)
f = prob.evaluate_l(xu, xl);
end

function [c, ce] = constraint_func(prob, xu, xl)
c = [];
ce = [];
end
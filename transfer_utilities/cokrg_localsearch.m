function [match_xl, lower_evalcount, trgdata, localsearch_fail, co_mdl] = cokrg_localsearch(xu, prob, cokrg_trg, lower_evalcount, maxFE)
%  This function does:
%  (1) builds cokrging model
%  (2) use cokrging to propose starting point (call: cokriging_decisionmaking)
%  (3) conduct local search from cokrging starting point (call: localsearch_withTransferredSolution)
%  Needs for record:
%  lower_evalcount -  how many evaluations done
%  trgdata - all solutions evaluated
% Needs to return:
% match_xl - found solution
% onbound - in (1) whether proposed solution is on HF boundary, if on
% boundary then switch to global search
%-------------------------------------------------


% number of evaluation
lower_evalcount = lower_evalcount + size(cokrg_trg.expensive_x, 1);

cheap_x = cokrg_trg.cheap_x;
cheap_f = cokrg_trg.cheap_f ;

[~, ia, ~] = unique(cheap_x, 'rows');    % check unique
cheap_x = cheap_x(ia, :);
cheap_f = cheap_f(ia, :);

% fprintf('cheap candidate go through distance check \n');
[cheap_x, cheap_f] = keepdistance_rawdata(cheap_x, cheap_f, prob.xl_bl, prob.xl_bu);

% eliminate repeated solutions;
x_trg{1,:} = cheap_x;                    % cheap
y_trg{1,:} = cheap_f;

expensive_x = cokrg_trg.expensive_x;
expensive_f = cokrg_trg.expensive_f;


[~, ia, ~]  = unique(expensive_x, 'rows');    % check unique
expensive_x = expensive_x(ia, :);
expensive_f = expensive_f(ia, :);

% fprintf('expensive candidate go through distance check \n');
[expensive_x, expensive_f] = keepdistance_rawdata(expensive_x, expensive_f, prob.xl_bl, prob.xl_bu);

x_trg{2,:} = expensive_x;
y_trg{2,:} = expensive_f;

% Fitting the cokring model
try
    co_mdl = oodacefit(x_trg, y_trg);
catch ME
    fprintf('--cokriging fails switch to replacement \n');
    co_mdl = replace_dace(expensive_x, expensive_f);
end
% [co_mdl, co_mdlparam] = cokrgmodel_building(x_trg, y_trg, prob);


cokrg_lb = min(x_trg{2,:}, [], 1); % use cheap because it expensive overlaps with cheap
cokrg_ub = max(x_trg{2,:}, [], 1);

bd_funh = @(x)onbound_withoutProblemBound(x, cokrg_lb, cokrg_ub, prob);
[cokrg_optxl, onbound_check] = cokriging_decisionmaking(co_mdl, prob, cokrg_lb, cokrg_ub, bd_funh);

if ~onbound_check
     [match_xl, lower_evalcount, trgdata] = cokrg_postprocess(cokrg_optxl, prob, xu, maxFE, lower_evalcount, expensive_x, expensive_f);
     localsearch_fail = false;
    return
else
    % as this step will later switch to ea search, trgdata is filled there

    % boundary process is: if hit boundary, re-evaluate
    % re-create starting point
    % if still hit boundary, check problem boundary
    % if on problem boundary, then accept it
    % else refuse the starting point onbound_check=false

    [cokrg_optfl, ~] = prob.evaluate_l(xu, cokrg_optxl);

    % re-evaluate x and f
    expensive_x = [expensive_x; cokrg_optxl];
    expensive_f = [expensive_f; cokrg_optfl];

    [expensive_x, expensive_f] = keepdistance_rawdata(expensive_x, expensive_f, prob.xl_bl, prob.xl_bu);

    x_trg{2,:} = expensive_x;
    y_trg{2,:} = expensive_f;

    % re--Fitting the cokring model
    try
        co_mdl = oodacefit(x_trg, y_trg);
    catch ME
        fprintf('--cokriging fails switch to replacement \n');
        co_mdl = replace_dace(expensive_x, expensive_f);
    end

    cokrg_lb = min(x_trg{2,:}, [], 1); % use cheap because it expensive overlaps with cheap
    cokrg_ub = max(x_trg{2,:}, [], 1);
    % repropose next point, if still on bound accept it
    bd_funh = @(x)onbound_considerProblemBound(x, cokrg_lb, cokrg_ub, prob);
    [cokrg_optxl, onbound_check] = cokriging_decisionmaking(co_mdl, prob, cokrg_lb, cokrg_ub, bd_funh);

    if onbound_check
        trgdata = 'place holder';
        match_xl = [];
        lower_evalcount = lower_evalcount + 1;    % one more evaluation 
        localsearch_fail = true;
    else
        lower_evalcount = lower_evalcount + 1;
        [match_xl, lower_evalcount, trgdata] = cokrg_postprocess(cokrg_optxl, prob, xu, maxFE, lower_evalcount, expensive_x, expensive_f);
        localsearch_fail = false;
        return

    end
end
end

function [match_xl, lower_evalcount, trgdata]= cokrg_postprocess(cokrg_optxl, prob, xu, maxFE, lower_evalcount, expensive_x, expensive_f)
% eliminate situations where proposed solution is not as good as
% existing one
cokrg_optfl = prob.evaluate_l(xu, cokrg_optxl);
if all(cokrg_optfl - expensive_f <= 0)
    fprintf('cokrg starting point accepted  \n');
else
    [~, id] = sort(expensive_f);
    cokrg_optxl = expensive_x(id(1), :);
end

% start local search from cokrging proposed point
% the one more evaluation (i.e. cokrg_optfl) is counted in local search 
[match_xl, lower_evalcount, x_visited, f_visited] = localsearch_withTransferredSolution(prob, xu, cokrg_optxl, maxFE, lower_evalcount);

%  TO BE NOTED: in local search output.funcCount is not the same number
%  as history.x, history.f Therefore  recorded FEs is not the same size
%  as collected  solutions in trgdata
trgdata = [expensive_x, expensive_f; x_visited, f_visited];
end

function [flag] = onbound_withoutProblemBound(x, lb, ub, prob)
% this function  check in normalized space
% whether x is on boundary

flag_up = false;
flag_down = false;

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


function [flag] = onbound_considerProblemBound(x, lb, ub, prob)
% this function  check in normalized space
% whether x is on boundary

flag_up = false;
flag_down = false;

x_norm = (x - prob.xl_bl) ./ (prob.xl_bu - prob.xl_bl);
lb_norm = (lb - prob.xl_bl) ./ (prob.xl_bu - prob.xl_bl);
ub_norm = (ub - prob.xl_bl) ./ (prob.xl_bu - prob.xl_bl);

prob_unorm = ones(1, prob.n_lvar);
prob_lnorm = zeros(1, prob.n_lvar);

check_ub = ub_norm - x_norm;
id_ub = check_ub < 1e-6;
if any(id_ub > 0) 
   check_probound = prob_unorm - x_norm;
   if any(check_probound < 1e-6)
       flag_up = false; % if the second time it is still on problem bound, then consider it not on bound
   else
       flag_up = true; 
   end
 
end

check_lb = x_norm - lb_norm;
id_lb = check_lb < 1e-6;
if any(id_lb > 0)
    check_probound = x_norm - prob_lnorm;
    if any(check_probound < 1e-6)
        flag_up = false; % if the second time it is also on problem bound, then consider it not on bound
    else
        flag_up = true;
    end
end

flag = flag_up|flag_down;

end 
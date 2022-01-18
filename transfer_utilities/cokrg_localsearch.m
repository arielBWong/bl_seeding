function [match_xl, lower_evalcount, trgdata, onbound] = cokrg_localsearch(xu, prob, cokrg_trg, lower_evalcount, maxFE)
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
co_mdl = oodacefit(x_trg, y_trg);

cokrg_lb = min(x_trg{2,:}, [], 1); % use cheap because it expensive overlaps with cheap
cokrg_ub = max(x_trg{2,:}, [], 1);

bd_funh = @(x)onbound(x, cokrg_lb, cokrg_ub, prob);
[cokrg_optxl, onbound_check] = cokriging_decisionmaking(co_mdl, prob, cokrg_lb, cokrg_ub, bd_funh);

if ~onbound_check
     [match_xl, lower_evalcount, trgdata] = cokrg_postprocess(cokrg_optxl, prob, xu, maxFE, lower_evalcount, expensive_x, expensive_f);
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
    co_mdl = oodacefit(x_trg, y_trg);

    cokrg_lb = min(x_trg{2,:}, [], 1); % use cheap because it expensive overlaps with cheap
    cokrg_ub = max(x_trg{2,:}, [], 1);
    % repropose next point, if still on bound accept it
    bd_funh = @(x)onbound_considerProblemBound(x, cokrg_lb, cokrg_ub, prob);
    [cokrg_optxl, onbound_check] = cokriging_decisionmaking(co_mdl, prob, cokrg_lb, cokrg_ub, bd_funh);

    if onbound_check
        trgdata = 'place holder';
        match_xl = cokrg_optxl;
        lower_evalcount = lower_evalcount + 1;
    else
        [match_xl, lower_evalcount, trgdata] = cokrg_postprocess(cokrg_optxl, prob, xu, maxFE, lower_evalcount, expensive_x, expensive_f);
        return

    end
end
end

function [match_xl, lower_evalcount, trgdata]= cokrg_postprocess(cokrg_optxl, prob, xu, maxFE, lower_evalcount, expensive_x, expensive_f)
% eliminate situations where proposed solution is not as good as
% existing one
cokrg_optfl = prob.evaluate_l(xu, cokrg_optxl);
if all(cokrg_optfl - expensive_f <= 0)
    fprintf('cokrg starting point is better than all training  \n');
else
    [~, id] = sort(expensive_f);
    cokrg_optxl = expensive_x(id(1), :);
end

% start local search from cokrging proposed point
[match_xl, lower_evalcount, x_visited, f_visited] = localsearch_withTransferredSolution(prob, xu, cokrg_optxl, maxFE, lower_evalcount);

%  TO BE NOTED: in local search output.funcCount is not the same number
%  as history.x, history.f Therefore  recorded FEs is not the same size
%  as collected  solutions in trgdata
trgdata = [cokrg_trg.expensive_x, cokrg_trg.expensive_f; x_visited, f_visited];
end
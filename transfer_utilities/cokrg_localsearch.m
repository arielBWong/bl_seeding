function [match_xl, lower_eval, lowerlocal_record, trgdata, co_mdl, close_optxlcokring, cokrg_lb, cokrg_ub, expensive_x] = cokrg_localsearch(xu, prob, close_optxu, close_optxl, expensive_x, ...
    lower_trg, cheap_x, cheap_f, lower_eval, lowerlocal_record, global_xl, global_fl, close_id, maxFE, sigma)
% This function is used to use close xu and close xl to do local search 
% return matching xl
%

seed_fl = prob.evaluate_l(close_optxu, close_optxl);
cheap_x = [cheap_x; close_optxl];            % cheap
cheap_f = [cheap_f; seed_fl];

[~, ia, ~] = unique(cheap_x, 'rows');    % check unique
cheap_x = cheap_x(ia, :);
cheap_f = cheap_f(ia, :);
fprintf('cheap candidate go through distance check \n');
[cheap_x, cheap_f] = keepdistance_rawdata(cheap_x, cheap_f, prob.xl_bl, prob.xl_bu);

% eliminate repeated solutions;
x_trg{1,:} = cheap_x;                    % cheap
y_trg{1,:} = cheap_f;


% a bit cheating on short cut
xuu = repmat(xu, size(expensive_x, 1), 1);
[expensive_f,~] = prob.evaluate_l(xuu, expensive_x);

seed_newfl  = prob.evaluate_l(xu, close_optxl);
expensive_x = [expensive_x; global_xl; close_optxl];
expensive_f = [expensive_f; global_fl; seed_newfl];

% number of evaluation
lower_eval = lower_eval+size(expensive_x, 1);


[~, ia, ~]  = unique(expensive_x, 'rows');    % check unique
expensive_x = expensive_x(ia, :);
expensive_f = expensive_f(ia, :);
fprintf('expensive candidate go through distance check \n');
[expensive_x, expensive_f] = keepdistance_rawdata(expensive_x, expensive_f, prob.xl_bl, prob.xl_bu);

x_trg{2,:} = expensive_x;
y_trg{2,:} = expensive_f;

% Fitting the cokring model
co_mdl = oodacefit(x_trg, y_trg);

cokrg_lb = min(x_trg{2,:} ,[], 1); % use cheap because it expensive overlaps with cheap
cokrg_ub = max(x_trg{2,:} ,[], 1);
[close_optxlcokring] = cokriging_decisionmaking(co_mdl, prob, cokrg_lb, cokrg_ub, 'hx',expensive_x, 'hf',expensive_f, 'lx',cheap_x ,'lf', cheap_f,...
    'boundary_updatesize',   floor(size(expensive_x, 1)/2),  'xu', xu, 'extend', false);

%  on boundary solutions are suspecious
if onboundary_check(close_optxlcokring, cokrg_lb, cokrg_ub, prob.xl_bl, prob.xl_bu)
    expensive_x = [expensive_x; close_optxlcokring];
    close_optfl  = prob.evaluate_l(xu, close_optxlcokring);
    expensive_f = [expensive_f; close_optfl];
    
    
    [~, ia, ~]  = unique(expensive_x, 'rows');    % check unique
    expensive_x = expensive_x(ia, :);
    expensive_f = expensive_f(ia, :);
    fprintf('expensive candidate go through distance check again \n');
    [expensive_x, expensive_f] = keepdistance_rawdata(expensive_x, expensive_f, prob.xl_bl, prob.xl_bu);

    
    x_trg{2,:} = expensive_x;
    y_trg{2,:} = expensive_f;
    
    % Fitting the cokring model
    co_mdl = oodacefit(x_trg, y_trg);
    
    cokrg_lb = min(x_trg{2,:} ,[], 1); % use cheap because it expensive overlaps with cheap
    cokrg_ub = max(x_trg{2,:} ,[], 1);
    [close_optxlcokring] = cokriging_decisionmaking(co_mdl, prob, cokrg_lb, cokrg_ub, 'hx',expensive_x, 'hf',expensive_f, 'lx',cheap_x ,'lf', cheap_f,...
        'boundary_updatesize',   floor(size(expensive_x, 1)/2),  'xu', xu, 'extend', false);
       
end


fl_optxlcokring = prob.evaluate_l(xu, close_optxlcokring);
if all(fl_optxlcokring - expensive_f <= 0)
    fprintf('cokrg starting point is better than all training \n');
else
    [~, id] = sort(expensive_f);
    close_optxlcokring = expensive_x(id(1), :);
end

% if after all process close_optxlcokring is still on boundary
% expand HF samples, how to expand? 
if onboundary_check(close_optxlcokring, cokrg_lb, cokrg_ub, prob.xl_bl, prob.xl_bu) && ~onboundary_check(close_optxlcokring, prob.xl_bl, prob.xl_bu, prob.xl_bl, prob.xl_bu)
    expensive_x = [expensive_x; close_optxlcokring];
    expensive_f = [expensive_f; fl_optxlcokring];
    nc = 20;
    close_optxlcokring = expandingHF(expensive_x, expensive_f, cheap_x, cheap_f, prob, xu, sigma, nc, close_optxu);
end

%maxFE = 950;
[match_xl, lower_eval, lowerlocal_record, x_visited, f_visited] = localsearch_withTransferredSolution(prob, xu, close_optxlcokring, maxFE, ...
    'lower_eval', lower_eval, 'lowerlocal_record', lowerlocal_record, ...
    'close_optxl', close_optxl,'close_optxu', close_optxu,'close_optxlcokring', close_optxlcokring  , ...
    'cheap_x',cheap_x ,'cheap_f', cheap_f , 'expensive_x', expensive_x, 'expensive_f',expensive_f, ...
    'co_mdl', co_mdl, ...
    'close_archive', lower_trg{close_id});



trgdata = [expensive_x, expensive_f; x_visited, f_visited];
end
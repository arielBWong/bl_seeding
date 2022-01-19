function[match_xl, mdl, trgdata, lower_searchSwitchFlag] = llmatch_trueEvaluation(upper_xu, prob, xl_probe, varargin)
% method of searching for a match xl for xu.
% Problem(Prob) definition require certain formation for bilevel problems
% evaluation method should be of  form 'evaluation_l(xu, xl)'
%--------------------------------------------------------------------------
% distribute parameters
% parse input
%-------
p = inputParser;
addRequired(p, 'upper_xu');
addRequired(p, 'prob');
addParameter(p, 'visualization', false);
addParameter(p, 'seeding_only', false);  % conduct seeding controller
addParameter(p, 'archive_xu', []);  % upper level all evaluated xu
addParameter(p, 'archive_xl', []);  % lower level matching xl, might larger than xu due to asynsic update
addParameter(p, 'seeding_strategy', 1); % 1 nearest neighbour 2 cokrging starting
addParameter(p, 'threshold', 0.9); % 1 nearest neighbour 2 cokrging starting
parse(p, upper_xu, prob, xl_probe, varargin{:});
% -----------------------------------

xu = p.Results.upper_xu;
visualization = p.Results.visualization;
seeding_only = p.Results.seeding_only;
archive_xu = p.Results.archive_xu;
archive_xl = p.Results.archive_xl;
prob = p.Results.prob;
seeding_strategy = p.Results.seeding_strategy;
thr = p.Results.threshold;
%------------------------------------------------

lower_searchSwitchFlag = 0;
global lower_eval
global lower_trg
%-------------------------------------------------

cokrg_samplesize = 20;
mdl = 'position_holder';
trgdata = 'position_holder';
initmatrix = []; %

funh_obj = @(x)objective_func(prob, xu, x);
funh_con = @(x)constraint_func(prob, xu, x);


%-- lower level initialization
% include seeding xl
if ~isempty(archive_xu) && seeding_only % first generation on the upper level use global search
    % ideal_xl = prob.get_xlprime(xu);
    maxFE = 1000;    % maxFE for local search
    
    if seeding_strategy == 1  % closest landscape
        [~, idx]= pdist2(archive_xu, xu, 'euclidean', 'Smallest', 1);  % this xu is upper level new infill xu, not added into archive_xu
        close_optxl =  archive_xl(idx, :);
        [match_xl, ~, history, output] = lowerlevel_fmincon(close_optxl, maxFE, prob.xl_bl, prob.xl_bu, funh_obj, funh_con);
        
        lower_eval = lower_eval + output.funcCount;     
        trgdata = [history.x, history.fval];
        lower_searchSwitchFlag = 1;
        return;
    end  

    if seeding_strategy == 2    % test neighbour
        % this method use cokriging sample to determine a starting point
        % and local search
        % cokriging sample both consider in std and outside std
        [expensive_x, expensive_f, cheap_x, cheap_f, correlation] = ...
            cokrg_trainingSampleExtension(xu, prob, archive_xu, archive_xl,  lower_trg, cokrg_samplesize);
        
        [neighbour_optxu, neighbour_optxl, ~] = retrieve_neighbour(xu, lower_trg, archive_xu, archive_xl);
        neighbour_optfl = prob.evaluate_l(neighbour_optxu, neighbour_optxl);
        current_fl = prob.evaluate_l(xu, neighbour_optxl);
        
        cokrg_trg.expensive_x = [expensive_x; neighbour_optxl];
        cokrg_trg.expensive_f = [expensive_f; current_fl];
        cokrg_trg.cheap_x = [cheap_x; neighbour_optxl];
        cokrg_trg.cheap_f = [cheap_f; neighbour_optfl];
        
        if correlation > thr
            % expensive x went through distance check
            maxFE = maxFE - size(cokrg_trg.expensive_x, 1);
            % expensive_x is counted inside cokrg_localsearch
            [match_xl, lower_eval, trgdata, localsearch_fail, cokrg_mdl] = cokrg_localsearch(xu, prob, cokrg_trg, lower_eval, maxFE);   
            
            % if cokrging propose a point on boundary
            if localsearch_fail
                % lower_eval has counted the match_xl in cokrg_localsearch
                % method
                initmatrix = [cokrg_trg.expensive_x; match_xl];
               %  algorithm will continue to ea search part
            else
                lower_searchSwitchFlag = 1;
                % test purpose
                f_test = prob.evaluate_l(xu, match_xl);
                xl_prime = prob.get_xlprime(xu);
                fl_prime = prob.evaluate_l(xu, xl_prime);

%                 if abs(fl_prime - f_test) > 0.01 % || visualization
%                     % plot 
%                     [close_optxu, close_optxl, ~] = retrieve_neighbour(xu, lower_trg, archive_xu, archive_xl);
%                    
%                     plot2d_withCokring(xu, prob.xl_bl, prob.xl_bu, match_xl, close_optxl, prob, close_optxu, cokrg_mdl,...
%                         cokrg_trg.cheap_x, cokrg_trg.cheap_f, cokrg_trg.expensive_x, cokrg_trg.expensive_f)
% 
%                      % plot2d_withCokring(xu, lb, ub, match_xl, close_optxl, prob, close_optxu, cokrg_mdl,...
%                             % co_cheapx, co_cheapy, co_expensivex, co_expensivey)
%                 end
                %---- 
                return;
            end
        else
            initmatrix = cokrg_trg.expensive_x;
        end
    end
    
   if seeding_strategy == 3    % test local search only 
        % this method use cokriging sample to determine a starting point
        % and local search
        % cokriging sample both consider in std and outside std
        [expensive_x, expensive_f, ~, ~, correlation] = ...
            cokrg_trainingSampleExtension(xu, prob, archive_xu, archive_xl,  lower_trg, cokrg_samplesize);
        
        maxFE = maxFE - size(expensive_x, 1);
        if correlation > thr
            % retrieve ID
            [~,~, close_id] = retrieve_neighbour(xu, lower_trg, archive_xu, archive_xl);
            close_optxl =  archive_xl(close_id, :);
            [match_xl, ~, history, output] = lowerlevel_fmincon(close_optxl, maxFE, prob.xl_bl, prob.xl_bu, funh_obj, funh_con);
            
            % record training data and evaluation usage
            trgdata = [expensive_x, expensive_f;  ...
                              history.x, history.fval];                          
             lower_eval = lower_eval + size(expensive_x, 1) + output.funcCount; 
              lower_searchSwitchFlag = 1;
             return;
        else
            initmatrix = expensive_x;
        end
   end


    
end


if visualization
    fighn            = figure(1);
end

% insert global search
% run to here means transfer fails or first generation

param.gen = 19;                    % design problem of gsolver, it should be 20
param.popsize = 50;
lb = prob.xl_bl;
ub = prob.xl_bu;
num_xvar = prob.n_lvar;
% initmatrix is set early in this method;

[best_x, best_f, best_c, archive_search, ~] = gsolver(funh_obj, num_xvar, lb, ub, initmatrix, funh_con, param, 'visualize', false);
lower_eval = lower_eval + (param.gen + 1) * param.popsize - size(initmatrix, 1); 

match_xl = best_x;
trgdata = [archive_search.sols(:, 2: prob.n_lvar + 1), archive_search.sols(:, prob.n_lvar + 2: end)];

end



function f = denormzscore(trainy, fnorm)
[~, y_mean, y_std] = zscore(trainy);
f = fnorm * y_std + y_mean;
end

function f = objective_func(prob, xu, xl)
f = prob.evaluate_l(xu, xl);
end

function [c, ce] = constraint_func(prob, xu, xl)
c = [];
ce = [];
end

function hitoptimal = hitoptimal_check(xu, xl, prob)
% function to check whehter x is in the vicinity of optimal
% 

ideal_xl = prob.get_xlprime(xu);

ideal_xl = (ideal_xl - prob.xl_bl)./(prob.xl_bu - prob.xl_bl);
xl = (xl - prob.xl_bl) ./ (prob.xl_bu - prob.xl_bl);

d = sqrt(sum((ideal_xl - xl).^2, 2));

if d < 0.01
    hitoptimal = true;
else
    hitoptimal = false;
end

end

function hitoptimal = hitoptimal_check2(xu, xl, prob)
% function to check whehter x is in the vicinity of optimal
% 

ideal_xl = prob.get_xlprime(xu);

ideal_xl = (ideal_xl - prob.xl_bl)./(prob.xl_bu - prob.xl_bl);
xl = (xl - prob.xl_bl) ./ (prob.xl_bu - prob.xl_bl);

d = sqrt(sum((ideal_xl - xl).^2, 2));

if d < 0.1
    hitoptimal = true;
else
    hitoptimal = false;
end

end


function which_localoptimal = checkwhich_localoptimal(xu, xl, prob)
% when this function is called, we already know that xl is local optimal
% 
othercenters = mp_othercenters(prob, xu);

[~, I] = pdist2(othercenters, xl, 'euclidean', 'Smallest', 1);
if I == 1 || I == 2
    which_localoptimal = 1; % forward and backward   
end

if I == 3 || I == 4
    which_localoptimal = 3; % left and right
end

end

function whetherBoundary = checkwhether_hitboundary(xl_startingpoint, xl_bl, xl_bu)
% when this function is called, we already know that xl is local optimal
whetherBoundary = 1;
if any( abs(xl_startingpoint - xl_bl) < 1e-3)
    whetherBoundary = 2;
    return;
end

if any(abs(xl_bu - xl_startingpoint) < 1e-3)
    whetherBoundary = 2;
    return
end

end


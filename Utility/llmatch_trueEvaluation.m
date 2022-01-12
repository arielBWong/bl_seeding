function[match_xl, mdl, trgdata] = llmatch_trueEvaluation(upper_xu, prob, xl_probe, varargin)
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
addRequired(p, 'xl_probe');
addParameter(p, 'visualization', false);
addParameter(p, 'seeding_only', false);  % conduct seeding controller
addParameter(p, 'archive_xu', []);  % upper level all evaluated xu
addParameter(p, 'archive_xl', []);  % lower level matching xl, might larger than xu due to asynsic update
addParameter(p, 'decision_making', true); % if this variable is set to false, means lower level is global search
addParameter(p, 'seeding_strategy', 1); % 1 nearest neighbour 2 cokrging starting
parse(p, upper_xu, prob, xl_probe, varargin{:});
% -----------------------------------

xu = p.Results.upper_xu;
visualization = p.Results.visualization;
seeding_only = p.Results.seeding_only;
lower_archive_xu =p.Results.lower_archive_xu;
archive_xu = p.Results.archive_xu;
init_size = p.Results.init_size;
xl_probe = p.Results.xl_probe;
prob = p.Results.prob;
archive_xl = p.Results.archive_xl;
decision_making = p.Results.decision_making;
seeding_strategy = p.Results.seeding_strategy;
%------------------------------------------------

global lower_eval
global lowerlocal_record
global lower_trg
global lowerlocal_recordWhichoptimal

%-------------------------------------------------

cokrg_samplesize = 20;
mdl = 'position_holder';
trgdata = 'position_holder';
initmatrix = []; %

funh_obj = @(x)objective_func(prob, xu, x);
funh_con = @(x)constraint_func(prob, xu, x);


%-- lower level initialization
global_halfsize = 10;
% include seeding xl
if ~isempty(archive_xu) && seeding_only % first generation on the upper level use global search
    % ideal_xl = prob.get_xlprime(xu);
    maxFE = 1000; % maxFE for local search
    
    if seeding_strategy == 1 % closest landscape
        [~, idx]= pdist2(archive_xu, xu, 'euclidean', 'Smallest', 1);  % this xu is upper level new infill xu, not added into archive_xu
        close_optxl =  archive_xl(idx, :);
        [match_xl, ~, history, output] = lowerlevel_fmincon(close_optxl, maxFE, prob.xl_bl, prob.xl_bu, funh_obj, funh_con);
        
        lower_eval = lower_eval + output.funcCount;     
        trgdata = [history.x, history.fval];
        return;
    end
    
       

    if seeding_strategy == 2 % test neighbour
        % this method investigates a different sampling strategy
        % assumption is that assume similar landscape, but moves, so it is
        % better sample solutions that has certain distance from neighbour
        % optimal
        [expensive_x, cheap_x, cheap_f, correlation, close_optxu, close_optxl, global_xl, global_fl, sigma, close_id] = ...
            cokrg_trainingSampleExtension(xu, archive_xu, lower_trg, prob, cokrg_samplesize, archive_xl);
        
        if correlation > 0.8
            % expensive x went through distance check
            [match_xl, lower_eval, ~, trgdata, co_mdl, cokrg_xl,  cokrg_lb, cokrg_ub, expensive_x] = cokrg_localsearch(xu, prob, close_optxu, close_optxl, expensive_x, ...
                lower_trg, cheap_x, cheap_f, lower_eval, [], global_xl, global_fl, close_id, maxFE, sigma);   
            
            nx = size(expensive_x, 1);
            expensive_f = prob.evaluate_l(repmat(xu, nx, 1), expensive_x);
            return;
        else
            initmatrix = expensive_x;
        end
    end
    
end


xudistance_decision = true;
visual_local = false;


if decision_making                                                         % if decision_making is set false, means that global search
    if ~isempty(lower_archive_xu) && ~isempty(archive_xu)                  % when archive_xu is passed in, means can do closeness check

        % only use xu distance to determine similar landscape
        if xudistance_decision
            
            [expensive_x, cheap_x, cheap_f, correlation, close_optxu, close_optxl, global_xl, global_fl] = cokrg_trainingExtraction(xu, archive_xu, lower_trg, prob, global_halfsize, archive_xl, global_halfflag);                
            % a repeat but have to have it step
            dist = pdist2(xu, archive_xu);                                 % this xu is upper level new infill xu, not added into archive_xu
            [~, idx] = sort(dist);
            
            for i = 1: length(idx)
                if size(lower_trg{idx(i)}, 1) > 500
                    close_id = idx(i);
                    break
                end
            end
            %
            if correlation > 0.6 % turn to local search
                [match_xl, lower_eval, lowerlocal_record, trgdata] = cokrg_localsearch(xu, prob, close_optxu, close_optxl, expensive_x, ...
                    lower_trg, cheap_x, cheap_f, lower_eval, lowerlocal_record, global_xl, global_fl, close_id);
 
                if visual_local
                    plot2d_withCokring(xu, prob.xl_bl, prob.xl_bu, match_xl, close_optxl, prob, close_optxu, co_mdl,  close_optxlcokring, ...
                        cheap_x, cheap_f, expensive_x, expensive_f);
                    
                end
                return;
            else
                fprintf('Roll back to global search \n');
                initmatrix = [expensive_x; global_xl];
            end
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
% initmatrix = [];

[best_x, best_f, best_c, archive_xu, ~] = gsolver(funh_obj, num_xvar, lb, ub, initmatrix, funh_con, param, 'visualize', false);
s = 1;

lower_eval = lower_eval + (param.gen + 1) * param.popsize - size(initmatrix, 1); % because before similarity, there is an adding evaluation step

match_xl = best_x;

if ~hitoptimal_check(xu, match_xl, prob)
    lowerlocal_record = [lowerlocal_record; 0]; % local search fail  
    lowerlocal_recordWhichoptimal = [lowerlocal_recordWhichoptimal; checkwhich_localoptimal(xu, match_xl, prob)];
    %  1 forward and backward   
    %  3 left and right
     
else
    lowerlocal_record = [lowerlocal_record; 1];         % local search success
    lowerlocal_recordWhichoptimal = [lowerlocal_recordWhichoptimal; 0];
end
    

aftersearch_collectx = archive_xu.sols(:, 2: prob.n_lvar + 1);
aftersearch_collectf = archive_xu.sols(:, prob.n_lvar + 2: end);  % only considers one objective
% add ons
trgdata = [aftersearch_collectx, aftersearch_collectf];

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


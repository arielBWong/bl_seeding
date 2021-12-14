function[match_xl, mdl, trgdata] = llmatch_trueEvaluation(upper_xu, prob, init_size, xl_probe, varargin)
% method of searching for a match xl for xu.
% Problem(Prob) definition require certain formation for bilevel problems
% evaluation method should be of  form 'evaluation_l(xu, xl)'
%--------------------------------------------------------------------------
% distribute parameters
%
% parse input
%-------
p = inputParser;
addRequired(p, 'upper_xu');
addRequired(p, 'prob');
addRequired(p, 'init_size');
addRequired(p, 'xl_probe');
addParameter(p, 'visualization', false);
addParameter(p, 'lower_archive_xu', []);                           % probing archive_xu
addParameter(p, 'seeding_only', false);                            % conduct seeding controller
addParameter(p, 'restartn', 0);                                    % additional restart times
addParameter(p, 'archive_xu', []);                                 % upper level all evaluated xu
addParameter(p, 'lower_xl', []);                                   % lower level matching xl, might larger than xu due to asynsic update
addParameter(p, 'decision_making', true);                          % if this variable is set to false, means lower level is global search
addParameter(p, 'global_halfflag', false);
addParameter(p, 'seed', -1);
addParameter(p, 'seeding_strategy', 1); % 1 naive starting 2 cokrging starting
parse(p, upper_xu, prob, init_size, xl_probe, varargin{:});
% -----------------------------------

xu = p.Results.upper_xu;
visualization = p.Results.visualization;
seeding_only = p.Results.seeding_only;
rn = p.Results.restartn;
lower_archive_xu =p.Results.lower_archive_xu;
archive_xu = p.Results.archive_xu;
init_size = p.Results.init_size;
xl_probe = p.Results.xl_probe;
prob = p.Results.prob;
lower_xl = p.Results.lower_xl;
decision_making = p.Results.decision_making;
seed = p.Results.seed;
global_halfflag = p.Results.global_halfflag;
seeding_strategy = p.Results.seeding_strategy;
%------------------------------------------------

global lower_eval
global lowerlocal_record
global lowerlocalsuccess_record2
global lower_trg
global lowerlocal_recordWhichoptimal
global lowerlocal_recordHitboundary;
%-------------------------------------------------


mdl = 'position_holder';
trgdata = 'position_holder';
initmatrix = []; %

funh_obj = @(x)objective_func(prob, xu, x);
funh_con = @(x)constraint_func(prob, xu, x);


%-- lower level initialization
% if ~isempty(lower_archive_xu)                                      % for lower level archive_xu
% train_xl = xl_probe;

global_halfsize = 10;
% include seeding xl
if ~isempty(archive_xu) && seeding_only                              % first generation on the upper level use global search
    ideal_xl      = prob.get_xlprime(xu);
    [ideal_fl, ~] = prob.evaluate_l(xu, ideal_xl);
      
    if rn>0 
        maxFE = 500; else
        maxFE = 1000;
    end
     
    if seeding_strategy == 1                                         % closest landscape
        dist = pdist2(xu , archive_xu);                              % this xu is upper level new infill xu, not added into archive_xu
        [~, idx] = sort(dist);
        % close_id = idx(1);
        for i = 1: length(idx)
            if size(lower_trg{idx(i)}, 1) > 500
                close_id = idx(i);
                break
            end
        end

        close_optxu = archive_xu(close_id, :);
        close_optxl =  lower_xl(close_id, :);
        [match_xl,~, history, output] = lowerlevel_fmincon(close_optxl, maxFE, prob.xl_bl, prob.xl_bu, funh_obj, funh_con);
        
        lower_eval = lower_eval + output.funcCount;
        [match_fl, match_cl] = prob.evaluate_l(xu, match_xl);         % additional lazy step, can be extracted from local search results        
        trgdata = [history.x, history.fval];
        return;
        
    end
    
    
    if seeding_strategy == 2 % cokrging   
        [expensive_x, cheap_x, cheap_f, correlation, close_optxu, close_optxl, global_xl, global_fl, sigma] = cokrg_trainingExtraction(xu, archive_xu, lower_trg, prob, global_halfsize, lower_xl, global_halfflag);
        % A repeat but have to have it step
        dist = pdist2(xu, archive_xu);                                 % this xu is upper level new infill xu, not added into archive_xu
        [~, idx] = sort(dist);
        
%         for i = 1: length(idx)
%             if size(lower_trg{idx(i)}, 1) > 500
%                 close_id = i;
%                 break
%             end
%         end
        close_id = idx(1);
        maxFE = maxFE - size(expensive_x, 1) + 1; 
        
        % expensive x went through distance check
        [match_xl, lower_eval, ~, trgdata, co_mdl, cokrg_xl,  cokrg_lb, cokrg_ub, expensive_x] = cokrg_localsearch(xu, prob, close_optxu, close_optxl, expensive_x, ...
            lower_trg, cheap_x, cheap_f, lower_eval, [], global_xl, global_fl, close_id, maxFE, sigma);
        
        if  ~hitoptimal_check2(xu, match_xl, prob)
            nex = size(expensive_x, 1);
            expensive_f = prob.evaluate_l(repmat(xu, nex, 1), expensive_x);
            plot2d_withCokring(xu, prob.xl_bl, prob.xl_bu, match_xl, close_optxl, prob, close_optxu, co_mdl, cokrg_xl,...
                            cheap_x,cheap_f, expensive_x, expensive_f, []);
        end
        [match_fl, match_cl] = prob.evaluate_l(xu, match_xl);                  % additional lazy step, can be extracted from local search results
    end    
    
    if seeding_strategy == 3                                                   % random
        rand_xl  = lhsdesign(1, prob.n_lvar, 'criterion','maximin','iterations',1000);
        rand_xl  = repmat(prob.xl_bl, 1, 1) ...
            + repmat((prob.xl_bu - prob.xl_bl), 1, 1) .* rand_xl;
        
        [match_xl,~, history, output] = lowerlevel_fmincon(rand_xl, maxFE, prob.xl_bl, prob.xl_bu, funh_obj, funh_con);        
        lower_eval = lower_eval + output.funcCount;
        [match_fl, match_cl] = prob.evaluate_l(xu, match_xl);                  % additional lazy step, can be extracted from local search results
        trgdata = [history.x, history.fval];
    end
    
    if seeding_strategy == 4
        % this method investigates a different sampling strategy
        % assumption is that assume similar landscape, but moves, so it is
        % better sample solutions that has certain distance from neighbour
        % optimal
        [expensive_x, cheap_x, cheap_f, correlation, close_optxu, close_optxl, global_xl, global_fl, sigma, close_id] = cokrg_trainingSampleExtension(xu, archive_xu, lower_trg, prob, global_halfsize, lower_xl);
        
        if correlation > 0.8
            % expensive x went through distance check
            [match_xl, lower_eval, ~, trgdata, co_mdl, cokrg_xl,  cokrg_lb, cokrg_ub, expensive_x] = cokrg_localsearch(xu, prob, close_optxu, close_optxl, expensive_x, ...
                lower_trg, cheap_x, cheap_f, lower_eval, [], global_xl, global_fl, close_id, maxFE, sigma);   
            
            nx = size(expensive_x, 1);
            expensive_f = prob.evaluate_l(repmat(xu, nx, 1), expensive_x);
            
            if ~hitoptimal_check2(xu, match_xl, prob)
                plot2d_withCokring(xu, prob.xl_bl, prob.xl_bu, match_xl, close_optxl, prob, close_optxu, co_mdl,  cokrg_xl, ...
                    cheap_x, cheap_f, expensive_x, expensive_f);
            end

            return;
        else
            initmatrix = expensive_x;
        end

        
    end
    
    
    % random restart other local seasrch points
    if rn>0 % number of restart       
        % random_restart
        two_checks = [];
        
        % first test the quality of previous transfer and local search
        % if abs(ideal_fl - match_fl) > 1   % this match_fl is from previous local search
        if ~hitoptimal_check(xu, match_xl, prob)
            two_checks(1) = 0;              % fail to identify global optimal
        else
            two_checks(1) = 1;              % success to identify global optimal
        end
        
        restart_xl  = lhsdesign(rn, prob.n_lvar, 'criterion','maximin','iterations',1000);
        restart_xl  = repmat(prob.xl_bl, rn, 1) ...
            + repmat((prob.xl_bu - prob.xl_bl), rn, 1) .* restart_xl;
        
        additional_searchxl = [];
        additional_searchfl = [];
        additional_searchcl = [];
        
        for i = 1:rn
            [match_xlrestart,~,history, output] = lowerlevel_fmincon(restart_xl, maxFE, prob.xl_bl, prob.xl_bu,funh_obj, funh_con);
            lower_eval = lower_eval + output.funcCount;
            
            [match_flrestart, match_clrestart] = prob.evaluate_l(xu, match_xlrestart);
            additional_searchxl = [additional_searchxl; match_xlrestart];
            additional_searchfl = [additional_searchfl; match_flrestart];
            additional_searchcl = [additional_searchcl; match_clrestart];
            
            trgdata = [trgdata;  history.x, history.fval];
        end
        
        % second test of the quality of random start and local search
        % if abs(ideal_fl - match_flrestart) > 1  % this match_fl is from previous local search
        if ~hitoptimal_check(xu, match_xl, prob)
            two_checks(2) = 0;           % fail to identify global optimal
        else
            two_checks(2) = 1;            % success to identify global optimal
        end
        lowerlocalsuccess_record2 = [lowerlocalsuccess_record2; two_checks];
        
        % select with previous search results with seeding
        additional_searchxl  = [additional_searchxl; match_xl];
        additional_searchfl  = [additional_searchfl; match_fl];
        additional_searchcl  = [additional_searchcl; match_cl];
        
        [match_xl, match_fl, ~, ~] =  localsolver_startselection(additional_searchxl, additional_searchfl, additional_searchcl);
                
    end
    
    % the above is local search invoked, so record number of local search
    % invoked for each xu
    % and success or not
    % if abs(ideal_fl - match_fl) > 1

%     if ~hitoptimal_check(xu, match_xl, prob)
%         lowerlocal_record = [lowerlocal_record; 0]; % local search fail
%         lowerlocal_recordWhichoptimal = [lowerlocal_recordWhichoptimal; checkwhich_localoptimal(xu, match_xl, prob)];
%          %  1 forward and backward   
%          %  3 left and right
%          if seeding_strategy == 2  % cokrg
%             lowerlocal_recordHitboundary = [lowerlocal_recordHitboundary; checkwhether_hitboundary(cokrg_xl,  cokrg_lb, cokrg_ub)];
%             
%             nex = size(expensive_x, 1);
%             expensive_f = prob.evaluate_l(repmat(xu, nex, 1), expensive_x);
%             
%             plot2d_withCokring(xu, prob.xl_bl, prob.xl_bu, match_xl, close_optxl, prob, close_optxu, co_mdl,  cokrg_xl, ...
%                         cheap_x, cheap_f, expensive_x, expensive_f);
%          
%             
%          end
%          
%          if seeding_strategy == 1  % neighbour optimal
%              % plot2d_withNeighbourOptimal(xu, prob.xl_bl, prob.xl_bu, match_xl, close_optxl, prob, close_optxu);
%             
%          end
%                            
%     else
%         lowerlocal_record = [lowerlocal_record; 1]; % local search success
%         lowerlocal_recordWhichoptimal = [lowerlocal_recordWhichoptimal; 0];
%         if seeding_strategy == 2
%             lowerlocal_recordHitboundary = [lowerlocal_recordHitboundary; 0];
%          end
%     end

    
    
end

% evaluate/get training fl from xu_init and train_xl
% compatible with non-constriant problem
% [train_fl, train_fc] = prob.evaluate_l(xu_init, train_xl);
% lower_archive_xu =[lower_archive_xu; train_fl'];
% lower_eval = lower_eval + size(train_fl, 1);

% decide whether skip infill steps
% correlation_decision = false;
% correctionFunc_decision = false;

xudistance_decision = true;
visual_local = false;


if decision_making                                                  % if decision_making is set false, means that global search
    if ~isempty(lower_archive_xu) && ~isempty(archive_xu)                 % when archive_xu is passed in, means can do closeness check

        % only use xu distance to determine similar landscape
        if xudistance_decision
            
            [expensive_x, cheap_x, cheap_f, correlation, close_optxu, close_optxl, global_xl, global_fl] = cokrg_trainingExtraction(xu, archive_xu, lower_trg, prob, global_halfsize, lower_xl, global_halfflag);                
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

[best_x, best_f, best_c, archive_xu ,~] = gsolver(funh_obj, num_xvar, lb, ub, initmatrix, funh_con, param, 'visualize', false);
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


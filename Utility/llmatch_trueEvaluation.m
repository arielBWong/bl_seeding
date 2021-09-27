function[match_xl] = llmatch_trueEvaluation(upper_xu, prob, init_size, xl_probe, varargin)
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
addRequired(p, 'init_size');
addRequired(p, 'xl_probe');
addParameter(p, 'visualization', false);
addParameter(p, 'seed_xl', []);
addParameter(p, 'lower_archive', []);                   %  probing archive
addParameter(p, 'seeding_only', false);               %  conduct seeding controller
addParameter(p, 'restartn', 0);                             %  additional restart times
addParameter(p, 'archive', []);                             %   upper level all evaluated xu
addParameter(p, 'lower_xl', []);                             %  lower level matching xl, might larger than xu due to asynsic update
parse(p, upper_xu, prob, init_size, xl_probe, varargin{:});

xu = p.Results.upper_xu;
seed_xl = p.Results.seed_xl;
visualization = p.Results.visualization;
seeding_only = p.Results.seeding_only;
rn = p.Results.restartn;
lower_archive =p.Results.lower_archive;
archive = p.Results.archive;
init_size = p.Results.init_size;
xl_probe = p.Results.xl_probe;
prob = p.Results.prob;
lower_xl = p.Results.lower_xl;
%------------------------------------------------

l_nvar                  = prob.n_lvar;
upper_bound      = prob.xl_bu;
lower_bound      = prob.xl_bl;
xu_init                 = repmat(xu, init_size, 1);

%-------------------------------------------------
opts = optimset('fmincon');
funh_obj = @(x)objective_func(prob, xu, x);
funh_con = @(x)constraint_func(prob, xu, x);


%-- lower level initialization
% if ~isempty(lower_archive)                                 % for lower level archive
    train_xl        = xl_probe;
% else
%     train_xl        = lhsdesign(init_size,l_nvar,'criterion','maximin','iterations',1000);
%     train_xl        = repmat(lower_bound, init_size, 1) ...
                                + repmat((upper_bound - lower_bound), init_size, 1) .* train_xl;
% end

% include seeding xl
if ~isempty(seed_xl)  && seeding_only                                 %  only local search with seeding
    if rn>0
        maxFE = 100; else
        maxFE = 200;
    end
    opts.Display = 'off';
    opts.MaxFunctionEvaluations = maxFE;
    [match_xl, ~, ~, ~] = fmincon(funh_obj, seed_xl, [], [],[], [],  ...
        prob.xl_bl, prob.xl_bu, funh_con,opts);
    
    [match_fl, match_cl] = prob.evaluate_l(xu, match_xl); % additional lazy step, can be extracted from local search results
    
    % random restart other local seasrch points
    if rn>0 % number of restart
        restart_xl  = lhsdesign(rn, prob.n_lvar, 'criterion','maximin','iterations',1000);
        restart_xl  = repmat(prob.xl_bl, rn, 1) ...
            + repmat((prob.xl_bu - prob.xl_bl), rn, 1) .* restart_xl;
        
        additional_searchxl = [];
        additional_searchfl = [];
        additional_searchcl = [];
        for i = 1:rn
            %[restart_fl, restart_cl] = prob.evaluate_l(xu, restart_xl(i, :));
            % [match_xlrestart, flag, num_eval] = ll_localsearch(restart_xl(i, :), restart_fl , restart_cl, true, xu, prob, maxFE);
            % n_fev = n_fev + num_eval + 1;
            
            opts.Display = 'off';
            opts.MaxFunctionEvaluations = maxFE;
            [match_xlrestart, ~, ~, ~] = fmincon(funh_obj, restart_xl, [], [],[], [],  ...
                prob.xl_bl, prob.xl_bu, funh_con,opts);
   
            [match_flrestart, match_clrestart] = prob.evaluate_l(xu, match_xlrestart);
            additional_searchxl = [additional_searchxl; match_xlrestart];
            additional_searchfl = [additional_searchfl; match_flrestart];
            additional_searchcl = [additional_searchcl; match_clrestart];
        end
        
        % select with previous search results with seeding
        additional_searchxl  = [additional_searchxl; match_xl];
        additional_searchfl   = [additional_searchfl; match_fl];
        additional_searchcl  = [additional_searchcl; match_cl];
        
        [match_xl, ~ , ~, ~] =  localsolver_startselection(additional_searchxl, additional_searchfl, additional_searchcl);
    end
    
    return
    
end

% evaluate/get training fl from xu_init and train_xl
% compatible with non-constriant problem
[train_fl, train_fc] = prob.evaluate_l(xu_init, train_xl);
lower_archive =[lower_archive; train_fl'];


% decide whether skip infill steps
if ~isempty(lower_archive) && ~isempty(archive)  % when archive is passed in, means can do closeness check
    dist  = pdist2(xu , archive);                              % this xu is upper level new infill xu, not added into archive
    [~, idx] = sort(dist);
    r = corr(train_fl, lower_archive(idx(1), :)');
    
    
    if r>0.95 % skip infill
        maxFE = 180;
        seed_xl = lower_xl(idx(1), :);
        [seed_fl, seed_fc] = prob.evaluate_l(xu, seed_xl);
        
        opts.Display = 'off';
        opts.MaxFunctionEvaluations = maxFE;
        [match_xl, ~, ~, ~] = fmincon(funh_obj, seed_xl, [], [],[], [],  ...
            prob.xl_bl, prob.xl_bu, funh_con,opts);        
        
        %  n_fev = num_eval + 1;
        [match_fl, match_cl]          = prob.evaluate_l(xu, match_xl);  % additional lazy step, can be extracted from local search results
        additional_searchxl           = [seed_xl; match_xl]; % variable name from copy paste
        additional_searchfl           = [seed_fl; match_fl];
        additional_searchcl           = [seed_fc; match_cl];
        
        [match_xl, ~ , flag, ~] =  localsolver_startselection(additional_searchxl, additional_searchfl, additional_searchcl);
       %  n_fev = n_fev + size(train_xl, 1);
        return
    end
end

if visualization
    fighn            = figure(1);
end


% insert global search
param.gen = 8;
param.popsize =  20;
lb = prob.xl_bl;
ub = prob.xl_bu;
num_xvar = prob.n_lvar;
initmatrix = train_xl;

[best_x, best_f, best_c, ~] = gsolver(funh_obj, num_xvar, lb, ub, initmatrix, funh_con, param);
s = 1;


nolocalsearch = false;
if nolocalsearch
    match_xl = best_x;
    n_fev      = size(arc_xl, 1);
    flag        = s;
else
    if size(train_fl, 2)> 1
        error('local search does not apply to MO');
    end
    maxFE = 20;
    % local search on true evaluation
  
    opts.Display = 'off';
    opts.MaxFunctionEvaluations = maxFE;
    
    [match_xl, ~, ~, ~] = fmincon(funh_obj, best_x, [], [],[], [],  ...
        prob.xl_bl, prob.xl_bu, funh_con,opts);
    
    % n_global                   = size(train_xl, 1);
    %  n_fev                      = n_global +num_eval;       % one in a population is evaluated
    
    % --- avoid sqp overshooting problem
    [match_fl, match_cl]         = prob.evaluate_l(xu, match_xl); % lazy step, no FE should be counted here
    additional_searchxl           = [best_x; match_xl]; % variable name from copy paste
    additional_searchfl           = [best_f;  match_fl];
    additional_searchcl           = [best_c; match_cl];
    
    [match_xl, ~ , flag, ~] =  localsolver_startselection(additional_searchxl, additional_searchfl, additional_searchcl);
    
end

% save lower level
% llcmp = true;
llcmp = false;
if llcmp
    % only for SO
    if size(train_fl, 2) ==  1
        arc_xl                = [arc_xl; match_xl];
        [local_fl, local_fc]  = prob.evaluate_l(xu, match_xl);
        arc_fl                = [arc_fl; local_fl];
        arc_cl                = [arc_cl; local_fc];
    end
    savelower(prob,arc_xl, arc_fl, arc_cl, method, seed);
end
end



function savelower(prob, x, f, c, method, seed)
num = length(prob.xl_bl);
savepath = strcat(pwd, '\resultfolder_ll\');
n = exist(savepath);
if n ~= 7
    mkdir(savepath)
end


savepath = strcat(pwd, '\resultfolder_ll\', prob.name, '_', num2str(num) ,'_',method);
n = exist(savepath);
if n ~= 7
    mkdir(savepath)
end

% check and process for multiple objective
% save nd front only
if size(f, 2)> 1
    % extract nd front
    num_con = size(c, 2);
    if ~isempty(c) % constraint problems
        index_c = sum(c <= 0, 2) == num_con;
        if sum(index_c) ~=0  % exist feasible,
            feasible_y = f(index_c, :);
            feasible_x = x(index_c, :);
            feasible_c = c(index_c, :);
            
            nd_index = Paretoset(feasible_y);
            f_nd = feasible_y(nd_index, :);
            x_nd = feasible_x(nd_index, :);
            c_nd = feasible_c(nd_index, :);
            
        else
            f_nd =[];
            x_nd = [];
            c_nd = [];
        end
    else % non constraint upper problem/ nd always exists
        nd_index = Paretoset(f);
        f_nd = f(nd_index, :);
        x_nd = x(nd_index, :);
        c_nd = [];
    end
    %---------------
    % save nd front
    savename_xu = strcat(savepath, '\xl_', num2str(seed),'.csv');
    savename_fu = strcat(savepath, '\fl_', num2str(seed),'.csv');
    savename_fc = strcat(savepath, '\cl_', num2str(seed),'.csv');
    
    csvwrite(savename_xu, x_nd);
    csvwrite(savename_fu, f_nd);
    csvwrite(savename_fc, c_nd);
    
    return
end

savename_xu = strcat(savepath, '\xl_', num2str(seed),'.csv');
savename_fu = strcat(savepath, '\fl_', num2str(seed),'.csv');
savename_fc = strcat(savepath, '\cl_', num2str(seed),'.csv');

csvwrite(savename_xu, x);
csvwrite(savename_fu, f);
csvwrite(savename_fc, c);
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


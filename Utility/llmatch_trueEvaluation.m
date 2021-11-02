function[match_xl, mdl] = llmatch_trueEvaluation(upper_xu, prob, init_size, xl_probe, varargin)
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
addParameter(p, 'lower_archive', []);                           % probing archive
addParameter(p, 'seeding_only', false);                         % conduct seeding controller
addParameter(p, 'restartn', 0);                                 % additional restart times
addParameter(p, 'archive', []);                                 % upper level all evaluated xu
addParameter(p, 'lower_xl', []);                                % lower level matching xl, might larger than xu due to asynsic update
addParameter(p, 'decision_making', true);                       % if this variable is set to false, means lower level is global search 
parse(p, upper_xu, prob, init_size, xl_probe, varargin{:});
% -----------------------------------

xu = p.Results.upper_xu;
visualization = p.Results.visualization;
seeding_only = p.Results.seeding_only;
rn = p.Results.restartn;
lower_archive =p.Results.lower_archive;
archive = p.Results.archive;
init_size = p.Results.init_size;
xl_probe = p.Results.xl_probe;
prob = p.Results.prob;
lower_xl = p.Results.lower_xl;
decision_making = p.Results.decision_making;
%------------------------------------------------

global lower_eval
global lowerlocal_record
global lower_mdl

l_nvar           = prob.n_lvar;
upper_bound      = prob.xl_bu;
lower_bound      = prob.xl_bl;
xu_init          = repmat(xu, init_size, 1);

%-------------------------------------------------


mdl = 'position_holder';


funh_obj = @(x)objective_func(prob, xu, x);
funh_con = @(x)constraint_func(prob, xu, x);


%-- lower level initialization
% if ~isempty(lower_archive)                                    % for lower level archive
train_xl = xl_probe;

% include seeding xl
if ~isempty(archive)  && seeding_only                           % first generation on the upper level use global search                                 
    
    if rn>0
        maxFE = 500; else
        maxFE = 1000;
    end
    dist = pdist2(xu , archive);                               % this xu is upper level new infill xu, not added into archive
    [~, idx] = sort(dist);
    seed_xl =  lower_xl(idx(1), :);
    
    % opts.Display = 'off';
    % opts.MaxFunctionEvaluations = maxFE;
    % [match_xl, ~, ~, output] = fmincon(funh_obj, seed_xl, [], [],[], [],  ...
    %     prob.xl_bl, prob.xl_bu, funh_con,opts);
    
    [match_xl,~, history, output] = lowerlevel_fmincon(seed_xl, maxFE, prob.xl_bl, prob.xl_bu, funh_obj, funh_con);
    
    lower_eval = lower_eval + output.funcCount;
    [match_fl, match_cl] = prob.evaluate_l(xu, match_xl); % additional lazy step, can be extracted from local search results
    
    % Collect solutions from all evaluated solutions
    % local search might not have as many as counted
    aftersearch_collectx = history.x;
    aftersearch_collectf = history.fval;
    
    % random restart other local seasrch points
    if rn>0 % number of restart
        restart_xl  = lhsdesign(rn, prob.n_lvar, 'criterion','maximin','iterations',1000);
        restart_xl  = repmat(prob.xl_bl, rn, 1) ...
            + repmat((prob.xl_bu - prob.xl_bl), rn, 1) .* restart_xl;
        
        additional_searchxl = [];
        additional_searchfl = [];
        additional_searchcl = [];
        
        for i = 1:rn           
            % opts.Display = 'off';
            % opts.MaxFunctionEvaluations = maxFE;  
            % [match_xlrestart, ~, ~, output] = fmincon(funh_obj, restart_xl, [], [],[], [],  ...
            %    prob.xl_bl, prob.xl_bu, funh_con,opts);
            
            [match_xlrestart,~,history, output] = lowerlevel_fmincon(restart_xl, maxFE, prob.xl_bl, prob.xl_bu,funh_obj, funh_con);
            
            lower_eval = lower_eval + output.funcCount;
            
            [match_flrestart, match_clrestart] = prob.evaluate_l(xu, match_xlrestart);
            additional_searchxl = [additional_searchxl; match_xlrestart];
            additional_searchfl = [additional_searchfl; match_flrestart];
            additional_searchcl = [additional_searchcl; match_clrestart];
            
            
            aftersearch_collectx = [aftersearch_collectx; history.x];
            aftersearch_collectf = [aftersearch_collectf; history.fval];
        end
        
        % select with previous search results with seeding
        additional_searchxl  = [additional_searchxl; match_xl];
        additional_searchfl  = [additional_searchfl; match_fl];
        additional_searchcl  = [additional_searchcl; match_cl];
        
        [match_xl, match_fl , ~, ~] =  localsolver_startselection(additional_searchxl, additional_searchfl, additional_searchcl);
    end 
    
    % the above is local search invoked, so record number of local search
    % invoked for each xu
    % and success or not
    
    ideal_xl = prob.get_xlprime(xu);
    [ideal_fl, ~] = prob.evaluate_l(xu, ideal_xl);
    if abs(ideal_fl - match_fl) > 1 
        lowerlocal_record = [lowerlocal_record;0]; % local search fail
    else
        lowerlocal_record = [lowerlocal_record;1]; % local search success
    end
    
    return
    
end

% evaluate/get training fl from xu_init and train_xl
% compatible with non-constriant problem
[train_fl, train_fc] = prob.evaluate_l(xu_init, train_xl);
% lower_archive =[lower_archive; train_fl'];
lower_eval = lower_eval + size(train_fl, 1);

% decide whether skip infill steps
correlation_decision = false;
correctionFunc_decision = true;
if decision_making                                                  % if decision_making is set false, means that global search
    if ~isempty(lower_archive) && ~isempty(archive)                 % when archive is passed in, means can do closeness check
        dist  = pdist2(xu, archive);                                % this xu is upper level new infill xu, not added into archive
        [~, idx] = sort(dist);
        r = corr(train_fl, lower_archive(idx(1), :)');
        
        seed_xu = archive(idx(1), :);
        
        if r > 0.95 % skip infill
            maxFE = 950;
            if correlation_decision
                seed_xl = lower_xl(idx(1), :);              
            end
            
            if correctionFunc_decision
                % Predict global optimal using landscape surrogate +
                % delta_D
         
                mdl_close = lower_mdl{idx(1)};         
                seed_xl = correctionFunc_decisionmaking(train_xl, train_fl,lower_archive(idx(1), :)', mdl_close, prob);
      
                
                
            end
            
            [seed_fl, seed_fc] = prob.evaluate_l(xu, seed_xl);
            [match_xl, ~, history, output] = lowerlevel_fmincon(seed_xl, maxFE, prob.xl_bl, prob.xl_bu,funh_obj, funh_con);
            
            lower_eval = lower_eval + output.funcCount;
            
            [match_fl, match_cl]          = prob.evaluate_l(xu, match_xl);  % additional lazy step, can be extracted from local search results
            additional_searchxl           = [seed_xl; match_xl]; % variable name from copy paste
            additional_searchfl           = [seed_fl; match_fl];
            additional_searchcl           = [seed_fc; match_cl];
            
            [match_xl, match_fl, ~, flag, ~] =  localsolver_startselection(additional_searchxl, additional_searchfl, additional_searchcl);
            
            ideal_xl = prob.get_xlprime(xu);
            [ideal_fl, ~] = prob.evaluate_l(xu, ideal_xl);
            if abs(ideal_fl - match_fl) > 1
                lowerlocal_record = [lowerlocal_record; 0]; % local search fail
            else
                lowerlocal_record = [lowerlocal_record; 1]; % local search success
            end
            
            % Collect solutions from all evaluated solutions
            % local search might not have as many as counted
            aftersearch_collectx = history.x;
            aftersearch_collectf = history.fval;
            
            aftersearch_collectx = [aftersearch_collectx; train_xl];
            aftersearch_collectf = [aftersearch_collectf; train_fl];
            mdl = lower_surrogateCreate(aftersearch_collectx, aftersearch_collectf, prob.xl_bl, prob.xl_bu);
            
            return
        end
    end
end


if visualization
    fighn            = figure(1);
end


% insert global search
% run to here means transfer fails or first generation

param.gen = 19; % design problem of gsolver, it should be 20
param.popsize = 50;
lb = prob.xl_bl;
ub = prob.xl_bu;
num_xvar = prob.n_lvar;
initmatrix = train_xl;

[best_x, best_f, best_c, archive ,~] = gsolver(funh_obj, num_xvar, lb, ub, initmatrix, funh_con, param, 'visualize', false);
s = 1;

lower_eval = lower_eval + (param.gen + 1) * param.popsize - size(train_xl, 1); % because before similarity, there is an adding evaluation step

match_xl = best_x;
if decision_making  % only when starting point selection is based on decision making, landscape surrogate needs to be trained and saved
    aftersearch_collectx = archive.sols(:, 2: prob.n_lvar + 1);
    aftersearch_collectf = archive.sols(:, prob.n_lvar + 2: end);  % only considers one objective
    mdl = lower_surrogateCreate(aftersearch_collectx, aftersearch_collectf, prob.xl_bl, prob.xl_bu);
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

function [f, gradient] = objective_func_extended(prob, xu, xl)
f = prob.evaluate_l(xu, xl);
gradient = prob.lower_gradient(xu, xl);

end


function plot2dlower(xu, lb, ub, match_xl, seed_xl, prob, seed_xu)

fignh = figure(3);
nt                  = 100;

cp                  = zeros(nt, nt);
x1_tst              = linspace(lb(1), ub(1), nt);
x2_tst              = linspace(lb(2), ub(2), nt);
[msx1, msx2]        = meshgrid(x1_tst, x2_tst);
msx11 = msx1(:);
msx22 = msx2(:);
xl = [msx11, msx22];
xum = repmat(xu, nt*nt, 1);
fl = prob.evaluate_l(xum, xl);
fl = reshape(fl, [nt, nt]);

fl_seed = prob.evaluate_l(xu, seed_xl);
fl_match = prob.evaluate_l(xu, match_xl);

subplot(1,2,1);
surf(msx1, msx2, fl, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;
scatter3(seed_xl(1), seed_xl(2), fl_seed,  80, 'r', 'filled' ); hold on;
scatter3(match_xl( 1), match_xl( 2), fl_match,  80, 'g', 'filled' ); hold on;

subplot(1,2,2);
xum = repmat(seed_xu, nt*nt, 1);
fl = prob.evaluate_l(xum, xl);
fl = reshape(fl, [nt, nt]);
fl_seed = prob.evaluate_l(seed_xu, seed_xl);
surf(msx1, msx2, fl, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;
scatter3(seed_xl(1), seed_xl(2), fl_seed,  80, 'r', 'filled' ); hold on;


pause(1);
close(fignh);


end

function localsearch_cheat(prob, seed_xl)
% this cheating solution is 
end


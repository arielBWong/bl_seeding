function blmapping_trueEvaldemo(problem_str, seed, varargin)
%  This function apply kriging to bl search and bl mapping
%
% parse input
p = inputParser;
addRequired(p, 'problem_str');
addRequired(p, 'seed');
addParameter(p, 'use_seeding', false);
addParameter(p, 'seeding_strategy', 1);
addParameter(p, 'threshold', 0.8);
parse(p, problem_str, seed, varargin{:});

%--- assign input parameters
prob_str = p.Results.problem_str;
seed = p.Results.seed;
use_seeding = p.Results.use_seeding;
seeding_strategy = p.Results.seeding_strategy;
thr = p.Results.threshold;
%---

visualize = false;
if visualize
    fighn  = figure('Position', [100 100 800 800]);
    h1 = subplot(2, 2, 1);
    h2 = subplot(2, 2, 3);
    h3 = subplot(2, 2, 4);
end

rng(seed, 'twister');
global prob
prob = eval(prob_str);
 
existing = result_check(prob, seed, use_seeding,  seeding_strategy);
if existing
    return;
end

%------------------Process starts--------------------
% insert global search

inisize_l = 20;
xl_probe = lhsdesign(inisize_l, prob.n_lvar, 'criterion','maximin','iterations',1000);
xl_probe = repmat(prob.xl_bl, inisize_l, 1) ...
    + repmat((prob.xl_bu - prob.xl_bl), inisize_l, 1) .* xl_probe;

funh_external = @(pop)up_probrecord(pop);
funh_obj  =  @(x)up_objective_func(prob, x, use_seeding, seeding_strategy, thr);
funh_con  =  @(x)up_constraint_func();

param.gen = 9;
param.popsize = 50;
lb = prob.xu_bl;
ub = prob.xu_bu;
num_xvar = prob.n_uvar;
initmatrix = [];

global upper_xu   % upper archive 
global lower_xl   % lower archive
global lower_eval
global lower_mdl
global lower_trg
global lower_decisionSwitch
global g
g = 1;

upper_xu = [];
lower_xl = [];

lower_eval = 0;
lower_mdl = {};
lower_trg = {};
lower_decisionSwitch = zeros(param.gen, param.popsize);

[best_x, ~, ~, ~, ~] = gsolver(funh_obj, num_xvar, lb, ub, initmatrix, funh_con, param,  'externalfunction', funh_external,  'visualize', false);


% Re-evaluation post process
% retrieve the last population
xu_lastpop = upper_xu(end - param.popsize + 1:end, :);
xl_lastpop = lower_xl(end - param.popsize + 1:end, :);
fu_lastpop = prob.evaluate_u(xu_lastpop, xl_lastpop);    % lazy step

lowertrgdata_lastpop = lower_trg(end - param.popsize + 1 : end);
indx_reeval = zeros(param.popsize, 1); % re-evaluation index

selected_xu = xu_lastpop(1, :);
selected_xl = xl_lastpop(1, :);

extra_lowerEval = 0;

reparam.gen = 19;
reparam.popsize = 50;

for ip = 1 : param.popsize
    % check condition (1) solution that is not fully evaluated
    % (2) after reevaluation position changes (i.e. no more best solution in the population)
    if indx_reeval(1,1) == 0 && size(lowertrgdata_lastpop{1}, 1) < 500 
        indx_reeval(1) = 1;
        
        
        xl_reeval = re_evaluation(prob, xu_lastpop(1, :), xl_lastpop(1, :), reparam); % always reevaluate top one
        
        extra_lowerEval = extra_lowerEval + reparam.popsize * (reparam.gen + 1) - 1;
        fu_reeval = prob.evaluate_u(xu_lastpop(1, :), xl_reeval);
        fu_lastpop(1, :) = fu_reeval;   % replace existing value;
        xl_lastpop(1, :) = xl_reeval;
        [~, idx] = sort(fu_lastpop);

        xu_lastpop = xu_lastpop(idx, : );
        xl_lastpop = xl_lastpop(idx, :);
        fu_lastpop = fu_lastpop(idx, :);
        lowertrgdata_lastpop = lowertrgdata_lastpop(idx);
        indx_reeval = indx_reeval(idx, :);
    else
        selected_xu = xu_lastpop(1, :);
        selected_xl = xl_lastpop(1, :); 
        break;
    end
end

save_results(upper_xu, lower_xl, prob,  selected_xu, selected_xl, seed, use_seeding, lower_eval, extra_lowerEval, seeding_strategy, lower_decisionSwitch, thr);

end

function xl_reeval = re_evaluation(prob, xu, xl, param)
% 
initmatrix = xl;
funh_obj = @(x)prob.evaluate_l(xu, x);
funh_con = @(x)up_constraint_func(); % re-use same outcome function

num_xvar = prob.n_lvar;
[xl_reeval, ~, ~, ~, ~] = gsolver(funh_obj, num_xvar, prob.xl_bl, prob.xl_bu, initmatrix, funh_con, param);

                            
end


function out = up_probrecord(pop)
% re-do one more round of initialization
% 

global upper_xu
global lower_xl
global lower_mdl
global lower_trg


upper_xu = [upper_xu; pop.X];
lower_xl = [lower_xl; pop.A];
lower_mdl = [lower_mdl, pop.Mdl];
lower_trg = [lower_trg, pop.trgdata];

out = [];
end

function [output] =  up_objective_func(prob, xu, use_seeding, seeding_strategy, thr)

global upper_xu
global lower_xl
global g
global lower_decisionSwitch
% upper xu and lower level does not change at the same time,
% upper_xu changes in generation wise, lower_xl changes in each xu's
% evaluation step. they should eventually have the same size.
% fprintf('upper generation: %d \n', g);

m = size(xu, 1);
f = [];
xl = [];
lower_searchSwitchFlags = [];
mdls = {};
trgdatas = {};

vis= false;
for i = 1:m
    fprintf('gen %d, ind %d \n ', g, i);
    xui = xu(i, :);

    if g == 3 && i == 14
        a = 0;
    end
    
    [match_xl, mdl, trgdata, lower_searchSwitchFlag] = llmatch_trueEvaluation(xui, prob, ...
        'archive_xu', upper_xu, 'archive_xl', lower_xl,...
        'seeding_only', use_seeding,  'seeding_strategy', seeding_strategy, ...
        'visualization', vis, 'threshold', thr);
    
    vis =  false;
    fi = prob.evaluate_u(xui, match_xl);
    
    xl = [xl; match_xl];
    f = [f; fi];
    mdls{end+1} = mdl;
    trgdatas{end+1} = trgdata;
    lower_searchSwitchFlags = [lower_searchSwitchFlags, lower_searchSwitchFlag];
end

fprintf('\n');
lower_decisionSwitch(g, :) = lower_searchSwitchFlags;
g = g + 1;

output.f = f;
output.addon = xl;
output.mdl = mdls;
output.trgdata = trgdatas;
end

function c = up_constraint_func()
c = [];
end



function existing = result_check(prob, seed, use_seeding,  seeding_strategy)
% this method is used for  breakpoint reconnect
existing = false;
name = strcat('resultfolder_trueEval', num2str(prob.n_lvar));
resultfolder = fullfile(pwd, name );

if use_seeding
    foldername = strcat(prob.name, '_seeding_strategy_', num2str(seeding_strategy));
else    
    foldername = strcat(prob.name, '_baseline_ea');
end

resultfolder = fullfile(resultfolder, foldername);

filename = strcat('xu_seed_', num2str(seed), '.csv');
savename = fullfile(resultfolder, filename);

if isfile(savename)
    existing = true;
end

end

function  save_results(xu, xl, prob, selected_xu, selected_xl, seed, use_seeding,  lower_eval, extra_lowerEval, seeding_strategy, lower_decisionSwitch, thr)

[fu, cu] = prob.evaluate_u(xu, xl);   % lazy  step
[fl, cl] = prob.evaluate_l(xu, xl);

thr = thr*100;
name = strcat('resultfolder_trueEval', num2str(prob.n_lvar), '_thr_', num2str(thr));
resultfolder = fullfile(pwd, name );
n = exist(resultfolder);
if n ~= 7
    mkdir(resultfolder)
end

if use_seeding
    foldername = strcat(prob.name, '_seeding_strategy_', num2str(seeding_strategy));
else    
    foldername = strcat(prob.name, '_baseline_ea');
end

resultfolder = fullfile(resultfolder, foldername);
n = exist(resultfolder);
if n ~= 7
    mkdir(resultfolder)
end

filename = strcat('xu_seed_', num2str(seed), '.csv');
savename = fullfile(resultfolder, filename);
csvwrite(savename, xu);

filename = strcat('xl_seed_', num2str(seed), '.csv');
savename = fullfile(resultfolder, filename);
csvwrite(savename, xl);

filename = strcat('fu_seed_', num2str(seed), '.csv');
savename = fullfile(resultfolder, filename);
csvwrite(savename, fu);

filename = strcat('fl_seed_', num2str(seed), '.csv');
savename = fullfile(resultfolder, filename);
csvwrite(savename, fl);


filename = strcat('lowerlevelcount_seed_', num2str(seed), '.csv');
savename = fullfile(resultfolder, filename);
csvwrite(savename, [lower_eval, extra_lowerEval]);

filename = strcat('lowerlevelswitch_seed_',  num2str(seed), '.mat');
savename = fullfile(resultfolder, filename);
save(savename,  'lower_decisionSwitch');

if size(fu, 2) > 1
    % fu nd front
    ndf_indx = Paretoset(fu);
    ndf = fu(ndf_indx, :);
    
    filename = strcat('final_ndfu_seed_', num2str(seed), '.csv');
    savename = fullfile(resultfolder, filename);
    csvwrite(savename, ndf);
    
    
    % fu hv
    pf = prob.upper_pf(100);
    max_f = max(pf, [], 1);
    min_f = min(pf, [], 1);
    norm_fu = (fu - min_f)./(max_f - min_f);
    ref = ones(1, size(pf, 2)) * 1.1;
    
    hv = Hypervolume(norm_fu, ref);
    filename = strcat('final_ndhv_seed_', num2str(seed), '.csv');
    savename = fullfile(resultfolder, filename);
    csvwrite(savename, hv);
    
    % fu hv
    igd = mean(min(pdist2(pf, ndf), [], 2));
    filename = strcat('final_ndigd_seed_', num2str(seed), '.csv');
    savename = fullfile(resultfolder, filename);
    csvwrite(savename, igd);
else
    % [best_x, best_f, best_c, s, index] =  localsolver_startselection(xu, fu, cu);
    best_fu = prob.evaluate_u(selected_xu, selected_xl);
    best_fl = prob.evaluate_l(selected_xu, selected_xl);
    ulp = abs(prob.fu_prime - best_fu);
    llp = abs(prob.fl_prime - best_fl);
    filename = strcat('final_accuracy_seed_', num2str(seed), '.csv');
    savename = fullfile(resultfolder, filename);
    csvwrite(savename, [ulp, llp]);
end


end

function [] = plotBothLevel(xui, match_xl, prob)

fignh               = figure(3);
fignh.Position      = [50, 100, 1000, 400];
nt                  = 100;

x1_tst              = linspace(prob.xu_bl(1), prob.xu_bu(1), nt);
x2_tst              = linspace(prob.xu_bl(2), prob.xu_bu(2), nt);
[msx1, msx2]        = meshgrid(x1_tst, x2_tst);
msx11 = msx1(:);
msx22 = msx2(:);
xu  = [msx11, msx22];
xlprime = prob.get_xlprime(xu);
fu  = prob.evaluate_u(xu, xlprime);
fu  = reshape(fu, [nt, nt]);

subplot(1, 2, 1);
surf(msx1, msx2, fu, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;

F = prob.evaluate_u(xui, match_xl);
scatter3(xui(1), xui(2), F,  80, 'r', 'filled' ); hold on;

subplot(1, 2, 2);
x1_tst              = linspace(prob.xl_bl(1), prob.xl_bu(1), nt);
x2_tst              = linspace(prob.xl_bl(2), prob.xl_bu(2), nt);
[msx1, msx2]        = meshgrid(x1_tst, x2_tst);
msx11 = msx1(:);
msx22 = msx2(:);
xl  = [msx11, msx22];
xu =  repmat(xui, nt*nt, 1);
f = prob.evaluate_l(xu, xl);
f = reshape(f, [nt, nt]);

surf(msx1, msx2, f,  'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;
f = prob.evaluate_l(xui, match_xl);
scatter3(match_xl(1), match_xl(2), f,  80, 'r', 'filled' ); hold on;

pause(1);
end



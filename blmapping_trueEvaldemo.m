function blmapping_trueEvaldemo(problem_str, seed, varargin)
%  This function apply kriging to bl search and bl mapping
%
% parse input
p = inputParser;
addRequired(p,'problem_str');
addRequired(p,'seed');
addParameter(p,'use_seeding', false);
addParameter(p,'seeding_strategy', 1);
parse(p, problem_str, seed, varargin{:});

%--- assign input parameters
prob_str = p.Results.problem_str;
seed = p.Results.seed;
use_seeding = p.Results.use_seeding;
seeding_strategy = p.Results.seeding_strategy;
%---

visualize = false;
if visualize
    fighn  = figure('Position', [100 100 800 800]);
    h1     = subplot(2, 2, 1);
    h2     = subplot(2, 2, 3);
    h3     = subplot(2, 2, 4);
end

rng(seed, 'twister');
global prob
prob = eval(prob_str);


%------------------Process starts--------------------
% insert global search
inisize_l = 20;
xl_probe = lhsdesign(inisize_l, prob.n_lvar, 'criterion','maximin','iterations',1000);
xl_probe = repmat(prob.xl_bl, inisize_l, 1) ...
    + repmat((prob.xl_bu - prob.xl_bl), inisize_l, 1) .* xl_probe;

funh_external = @(pop)up_probrecord(pop,  xl_probe, prob);


funh_obj  =  @(x)up_objective_func(prob, x, use_seeding, seeding_strategy);
funh_con  =  @(x)up_constraint_func(prob, x);

param.gen = 9;
param.popsize = 50;
lb = prob.xu_bl;
ub = prob.xu_bu;
num_xvar = prob.n_uvar;
initmatrix = [];

global upper_xu  % upper archive 
global lower_xl   % lower archive
global lower_eval
global lower_mdl
global lower_trg
global g
global gpr_mdl
g = 1;

upper_xu = [];
lower_xl = [];

lower_eval = 0;
lower_mdl = {};
lower_trg ={};
gpr_mdl = [];


[best_x, ~, ~, a, ~] = gsolver(funh_obj, num_xvar, lb, ub, initmatrix, funh_con, param,  'externalfunction', funh_external,  'visualize', false);


save_results(upper_xu, lower_xl, prob,  seed, use_seeding, lower_eval, seeding_strategy);
end

function out = up_probrecord(pop, xl_probe, prob)
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


function [output] =  up_objective_func(prob, xu, use_seeding, seeding_strategy)

global upper_xu
global lower_xl
global g
% upper xu and lower level does not change at the same time,
% upper_xu changes in generation wise, lower_xl changes in each xu's
% evaluation step. they should eventually have the same size.
% fprintf('upper generation: %d \n', g);

m = size(xu, 1);
f = [];
xl = [];
mdls = {};
trgdatas = {};

vis= false;
for i = 1:m
    fprintf('gen %d, ind %d \n ', g, i);
    xui = xu(i, :);
    
    [match_xl, mdl, trgdata] = llmatch_trueEvaluation(xui, prob, ...
        'archive_xu', upper_xu, 'archive_xl', lower_xl,...
        'seeding_only', use_seeding,  'seeding_strategy', seeding_strategy, ...
        'visualization', vis);
    
    fi = prob.evaluate_u(xui, match_xl);
    
    xl = [xl; match_xl];
    f = [f; fi];
    mdls{end+1} = mdl;
    trgdatas{end+1} = trgdata;
end

fprintf('\n');
g = g + 1;

output.f = f;
output.addon = xl;
output.mdl = mdls;
output.trgdata = trgdatas;
end


function c = up_constraint_func(prob, xu)
c = [];
end




function  save_results(xu, xl, prob, seed, use_seeding,  lower_eval,seeding_strategy)

[fu, cu] = prob.evaluate_u(xu, xl);   % lazy  step
[fl, cl] = prob.evaluate_l(xu, xl);

name = strcat('resultfolder_trueEval', num2str(prob.n_lvar));
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
csvwrite(savename, lower_eval);


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
    [best_x, best_f, best_c, s, index] =  localsolver_startselection(xu, fu, cu);
    ulp = abs(prob.fu_prime - best_f);
    llp = abs(prob.fl_prime - fl(index) );
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



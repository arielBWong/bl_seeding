function blmapping_trueEvaldemo(problem_str, seed, varargin)
%  This function apply kriging to bl search and bl mapping
%

% parse input
p = inputParser;
addRequired(p,  'problem_str');
addRequired(p,  'seed');
addParameter(p, 'restart_num', 0);
addParameter(p, 'use_seeding', false);
parse(p, problem_str, seed, varargin{:});

%--- assign input parameters
prob_str = p.Results.problem_str;
seed = p.Results.seed;
use_seeding = p.Results.use_seeding;
rx = p.Results.restart_num;
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
inisize_l       = 20;
xl_probe        = lhsdesign(inisize_l, prob.n_lvar, 'criterion','maximin','iterations',1000);
xl_probe        = repmat(prob.xl_bl, inisize_l, 1) ...
    + repmat((prob.xl_bu - prob.xl_bl), inisize_l, 1) .* xl_probe;

funh_external = @(pop)up_probrecord(pop,  xl_probe, prob);


funh_obj  =  @(x)up_objective_func(prob, x, xl_probe, use_seeding);
funh_con  =  @(x)up_constraint_func(prob, x);

param.gen = 10;
param.popsize = 50;
lb = prob.xu_bl;
ub = prob.xu_bu;
num_xvar = prob.n_uvar;
initmatrix = [];

global upper_xu
global xu_probefl
global lower_xl;
global lower_eval;
upper_xu = [];
lower_xl = [];
xu_probefl = [];
lower_eval = 0;

% obj = VideoWriter('moving.avi');
% obj.Quality= 100;
% obj.FrameRate = 25;
% open(obj);

[best_x, ~, ~, ~] = gsolver(funh_obj, num_xvar, lb, ub, initmatrix, funh_con, param,  'externalfunction', funh_external,  'visualize', false);


%obj.close();
save_results(upper_xu, lower_xl, prob,  seed, use_seeding, rx, lower_eval);
end

function out = up_probrecord(pop, xl_probe, prob)
% re-do one more round of initialization
%

global xu_probefl
global upper_xu
global lower_xl

upper_xu = [upper_xu; pop.X];
lower_xl = [lower_xl; pop.A];
n = size(xl_probe, 1);
m = size(pop.X, 1);
for i =1:m
    xui = repmat(pop.X(i, :), n, 1);
    fl_probe = prob.evaluate_l(xui, xl_probe);
    xu_probefl = [xu_probefl; fl_probe'];
end
out= [];
end


function [output] =  up_objective_func(prob, xu, xl_probe, use_seeding)
global xu_probefl
global upper_xu
global lower_xl

% upper xu and lower level does not change at the same time, 
% upper_xu changes in generation wise, lower_xl changes in each xu's
% evaluation step. they should eventually have the same size.


m = size(xu, 1);
f = [];
xl = [];
for i = 1:m
    xui = xu(i, :);
  
    [match_xl] =  llmatch_trueEvaluation(xui, prob,  20, xl_probe, ...
        'lower_archive', xu_probefl, 'archive', upper_xu, 'lower_xl', lower_xl,...
        'seeding_only', use_seeding);
    fi = prob.evaluate_u(xui, match_xl);
    xl = [xl; match_xl];
    f = [f; fi];
end
output.f = f;
output.addon = xl; 

end

function c = up_constraint_func(prob, xu)
c = [];
end

function  save_results(xu, xl, prob, seed, use_seeding, rx, lower_eval)

[fu, cu] = prob.evaluate_u(xu, xl); % lazy  step
[fl, cl] = prob.evaluate_l(xu, xl);


resultfolder = fullfile(pwd, 'resultfolder_trueEval' );
n = exist(resultfolder);
if n ~= 7
    mkdir(resultfolder)
end

if use_seeding
    foldername = strcat(prob.name, '_seeding');
    foldername = strcat(prob.name, num2str(rx));
else
    % foldername = strcat(prob.name, '_global2');
    foldername = prob.name;
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

filename = strcat('cu_seed_', num2str(seed), '.csv');
savename = fullfile(resultfolder, filename);
csvwrite(savename, cu);

filename = strcat('cl_seed_', num2str(seed), '.csv');
savename = fullfile(resultfolder, filename);
csvwrite(savename, cl);

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





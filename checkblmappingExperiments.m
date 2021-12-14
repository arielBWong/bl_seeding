% This code is related to design variable  mapping 
% This code tests whether using two points (BLmapping & nearest neighbour can in general have the same starting point)
%

% (1) Initialize 50 points for each problem
% (2) Find their lower level problem xl
% (3) Create two points one from close neighbour one from BL mapping
% (4) Create distance measure, check the euclidean distance


clear;
close all;

problem_folder = strcat(pwd, '/TP3');
addpath(problem_folder);
problem_folder = strcat(pwd, '/globalsolver');
addpath(problem_folder);
problem_folder = strcat(pwd, '/ND_Sort');
addpath(problem_folder);
problem_folder = strcat(pwd, '/OODACE');
addpath(problem_folder);


problems = { 'smd1mp(1, 2, 1)' , 'smd2mp(1, 2, 1)',  'smd3mp(1, 2, 1)',  'smd4mp(1, 2, 1)', ....
   'smd5mp(1, 2, 1)' , 'smd6mp(1, 0, 2, 1)', 'smd7mp(1, 2, 1)',  'smd8mp(1, 2, 1)'};

problems = { 'smd1mp(1, 1, 1)' , 'smd2mp(1, 1, 1)',  'smd3mp(1, 1, 1)',  'smd4mp(1, 1, 1)', ....
   'smd5mp(1, 1, 1)' , 'smd6mp(1, 0, 1, 1)', 'smd7mp(1, 1, 1)',  'smd8mp(1, 1, 1)'};

prob = smd5mp(1, 2, 1);
% checkblmapping_experiments(prob, 1);

save2documents(prob, 1, 1);

for seed = 1 : 29    
    for i = 1 : 8
        fprintf('%d\n', i);
        prob = eval(problems{i});
        checkblmapping_experiments(prob, seed);
    end   
end

function checkblmapping_experiments(prob, seed)
rng(seed, 'twister');
% (1)

nu = 50;
xu = lhsdesign(nu, prob.n_uvar, 'criterion', 'maximin', 'iterations', 1000);
xu = repmat(prob.xu_bl, nu, 1) + repmat((prob.xu_bu - prob.xu_bl), nu, 1) .* xu;

% (2) 
xl = [];
param.gen = 9;
param.popsize = 50;


for i = 1:nu
    funh_obj = @(x)lower_objective(prob, xu(i, :), x);
    funh_con = @(x)lower_constraint(prob, xu(i, :), x);
    
    [best_x, best_f, best_c, archive ,~] = gsolver(funh_obj, prob.n_lvar, prob.xl_bl, prob.xl_bu, [], funh_con, param, 'visualize', false);
    xl = [xl; best_x];
end


mdls = {};
for i = 1 : prob.n_lvar
    mdls{i} = oodacefit(xu, xl(:, i));
end

%(3-1)
ns = 50;
xu_tst = lhsdesign(ns, prob.n_uvar, 'criterion', 'maximin', 'iterations', 1000);
xu_tst = repmat(prob.xu_bl, ns, 1) + repmat((prob.xu_bu - prob.xu_bl), ns, 1) .* xu_tst;

%(4-1) check prediction
bl_startingpoint = {};
for k = 1:ns    
    xltst = [];
    for i = 1:prob.n_lvar
        tmp = mdls{i}.predict(xu_tst(k, :));
        xltst = [xltst, tmp];
    end   
    bl_startingpoint{k} = xltst;
end

% (4-2) check neighbour
neighbour_startingpoint = {};
for k = 1:ns
    [~, I] = pdist2(xu, xu_tst(k, :),  'euclidean', 'Smallest', 1); % process one by one
    xu_neighbour = xu(I, :);
    xl_neighbour = xl(I, :);
    neighbour_startingpoint{k} = xl_neighbour;
end

% (4-3) check second neighbour
neighbour2_startingpoint = {};
for k = 1:ns
    [~, I] = pdist2(xu, xu_tst(k, :),  'euclidean', 'Smallest', 2); % process one by one
    xu_neighbour2 = xu(I(2), :);
    xl_neighbour2 = xl(I(2), :);
    neighbour2_startingpoint{k} = xl_neighbour2;
end

% (5) compare starting point
% check 
for k = 1:ns
    f_map = prob.evaluate_l(xu_tst(k, :), bl_startingpoint{k});
    
    f_neighbour = prob.evaluate_l(xu_tst(k, :),  neighbour_startingpoint{k});    
    f_neighbour2 = prob.evaluate_l(xu_tst(k, :),  neighbour2_startingpoint{k});
    
    [~, idx] = min([f_map, f_neighbour], [], 2);
    if idx == 1
        selectedx{k} = bl_startingpoint{k};
    else
        selectedx{k} = neighbour_startingpoint{k};
    end
     
    [~, idx] = min([f_neighbour, f_neighbour2], [], 2);
    if idx == 1
        comparex{k} = neighbour_startingpoint{k};
    else
        comparex{k} = neighbour2_startingpoint{k};
    end
    %---------------   
end


% 
% 
d1 = [];
d2 = [];
d3 = [];
for k = 1:ns
    xl_prime = prob.get_xlprime(xu_tst(k, :));
    d1 = [d1; pdist2(xl_prime, selectedx{k}, 'euclidean')]; 
    d2 = [d2; pdist2(xl_prime, neighbour_startingpoint{k}, 'euclidean')];
    d3 = [d3; pdist2(xl_prime, comparex{k}, 'euclidean')];
end

% conduct significant test cu
% 1 means  the former is smaller,  second is bigger/better,
% for distance comparison smaller better

[p1,h1,stats1] = ranksum(d1, d3,  'alpha', 0.05, 'tail', 'left');
[p2,h2,stats2] = ranksum(d3, d1,  'alpha', 0.05, 'tail', 'left');
if h1 == 1 && h2 == 0
    fprintf('mapping+ method is significant smaller/better distance to optimal \n');
elseif h2 == 1 && h1 == 0
    fprintf('mapping+ method is worse \n');
else
    fprintf('no difference \n');
end


end


function save2documents(prob, seed, sigtest)
% 
%
%
foldername = 'result_folder3';
resultfolder = fullfile(pwd, foldername);

problem_folder = prob.name;
problem_folder = strcat(foldername, '\\',problem_folder);

if ~exist(problem_folder, 'dir')
    mkdir(problem_folder);
end

%




end


function obj = lower_objective(prob, xu, xl)
obj = prob.evaluate_l(xu, xl);
end

function con = lower_constraint(prob,xu, xl)
con = [];
end

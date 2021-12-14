clear;
close all;

problem_folder = strcat(pwd,'/problems/TP3');
addpath(problem_folder);

problem_folder = strcat(pwd,'/Surrogate/Methods/Surrogate');
addpath(problem_folder);

problem_folder = strcat(pwd,'/globalsolver');
addpath(problem_folder);

problem_folder = strcat(pwd,'/ND_Sort');
addpath(problem_folder);

problem_folder = strcat(pwd,'/Utility');
addpath(problem_folder);

problem_folder = strcat(pwd,'/OODACE');
addpath(problem_folder);

problem_folder = strcat(pwd,'/transfer_utilities');
addpath(problem_folder);

problems = {'smd5mp(1, 1, 1)'};

% tic;
cokrg_test('smd7mp(1, 1, 1)',  11, 'decision_making', false, 'use_seeding', true, 'restart_num', 0, 'seeding_strategy', 2);
% toc;


function cokrg_test(problem_str, seed, varargin)
%  This function is used to test why cokrg has performance that contradict
%  to lower level success rate
%
% parse input
p = inputParser;
addRequired(p,  'problem_str');
addRequired(p,  'seed');
addParameter(p, 'decision_making', true);
addParameter(p, 'restart_num', 0);
addParameter(p, 'use_seeding', false);
addParameter(p, 'seeding_strategy', 1);
parse(p, problem_str, seed, varargin{:});

%--- assign input parameters
prob_str = p.Results.problem_str;
seed = p.Results.seed;
use_seeding = p.Results.use_seeding;
rn = p.Results.restart_num;
decision_making = p.Results.decision_making;
seeding_strategy = p.Results.seeding_strategy;
%---

% generate 50 xu and apply EA to create corresponding lower level
% solutions
inisize = 50;
xu_base = lhsdesign(inisize, prob.n_uvar, 'criterion', 'maximin', 'iterations', 1000);
xu_base = repmat(prob.xu_bl, inisize, 1) ...
                + repmat((prob.xu_bu - prob.xu_bl), inisize, 1) .* xu_base;
output = up_objective_func(prob, xu_base, xl_probe, use_seeding, rn, decision_making, seed, seeding_strategy);


% generate another 100 solutions and apply different scheme (neighbour optimal + local / cokrg optimal + local)
inisize_new = 100;
xu_new = lhsdesign(inisize_new, prob.n_uvar, 'criterion', 'maximin', 'iterations', 1000);
xu_new = repmat(prob.xu_bl, inisize_new, 1) ...
                + repmat((prob.xu_bu - prob.xu_bl), inisize_new, 1) .* xu_new;
xl_prob = [];
output = up_objective_func(prob, xu, xl_probe, use_seeding, rn, decision_making, seed, seeding_strategy);


% save accuracy and success rate for comparison



end



function [output] =  up_objective_func(prob, xu, xl_probe, use_seeding, rn, decision_making, seed, seeding_strategy)
global xu_probefl
global upper_xu
global lower_xl
global g
% upper xu and lower level does not change at the same time, 
% upper_xu changes in generation wise, lower_xl changes in each xu's
% evaluation step. they should eventually have the same size.
fprintf('upper generation: %d \n', g);


g = g + 1;
m = size(xu, 1);
f = [];
xl = [];
mdls = {};
trgdatas = {};

vis= false;
for i = 1:m
    fprintf('ind %d \n ', i);
    xui = xu(i, :);
    
    if g == 5 && i == 16
        vis = true;
    end
    [match_xl, mdl, trgdata] = llmatch_trueEvaluation(xui, prob,  20, xl_probe, ...
        'lower_archive', xu_probefl, 'archive', upper_xu, 'lower_xl', lower_xl,...
        'seeding_only', use_seeding, 'restartn', rn, 'decision_making', decision_making,...
        'seed', seed, 'global_half', false, ...
        'seeding_strategy', seeding_strategy, ...
         'visualization', vis);
    
    vis = false;
    fi = prob.evaluate_u(xui, match_xl);
    
    if fi < -10
        plotBothLevel(xui, match_xl, prob);
    end
    xl = [xl; match_xl];
    f  = [f; fi];
    mdls{end+1} = mdl;
    trgdatas{end+1} = trgdata;
    
end
fprintf('\n');

output.f = f;
output.addon = xl; 
output.mdl = mdls;
output.trgdata = trgdatas;

end
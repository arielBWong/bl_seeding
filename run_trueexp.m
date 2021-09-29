%%this script is to run ulego with multiple seeds
%
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

problems = {'smd5mp(1, 1, 1)' , 'smd7mp(1, 1, 1)',  'smd8mp(1, 1, 1)'};
% tic;
% blmapping_trueEvaldemo('smd7mp(1, 1, 1)',  2, 'use_seeding', false, 'restart_num', 1);
% toc;

% --------
seeds = 1:29;
ns = length(seeds);
np = length(problems);
%------

% create parameter for each method
paras  = struct([]);

% create parameter for the first method
for i = 1 : np
    for j = 1: ns
        strc_id = (i-1) * ns + j;
        paras(strc_id). problem_str = problems{i};
        paras(strc_id). seed = j;
        paras(strc_id). use_seeding = false;
        paras(strc_id). restart_num = 0;
    end
end

% create parameter for the second method
for i = 1 : np
    for j = 1: ns
        strc_id = (i-1) * ns + j + ns *np;
        paras(strc_id). problem_str = problems{i};
        paras(strc_id). seed = j;
        paras(strc_id). use_seeding = true;
        paras(strc_id). restart_num = 0;
    end
end

% create parameter for the third method
for i = 1 : np
    for j = 1: ns
        strc_id = (i-1) * ns + j +ns *np*2;
        paras(strc_id). problem_str = problems{i};
        paras(strc_id). seed = j;
        paras(strc_id). use_seeding = true;
        paras(strc_id). restart_num = 1;
    end
end
nrun = length(paras);

parfor i = 1:nrun
    blmapping_trueEvaldemo(paras(i).problem_str, paras(i).seed,...
                             'use_seeding', paras(i).use_seeding, 'restart_num', paras(i).restart_num);
                        
end


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

problem_folder = strcat(pwd,'/OODACE');
addpath(problem_folder);

problem_folder = strcat(pwd,'/transfer_utilities');
addpath(problem_folder);

problems = {'smd5mp(1, 1, 1)' , 'smd7mp(1, 1, 1)',  'smd8mp(1, 1, 1)', ...
    'smd1mp(1, 1, 1)' , 'smd2mp(1, 1, 1)',  'smd3mp(1, 1, 1)', ...
    'smd4mp(1, 1, 1)', 'smd6mp(1, 0, 1, 1)', };


% problems = {'smd5mp(1, 2, 1)' , 'smd7mp(1, 2, 1)',  'smd8mp(1, 2, 1)', ...
%     'smd1mp(1, 2 , 1)' , 'smd2mp(1, 2, 1)',  'smd3mp(1, 2, 1)', ...
%     'smd4mp(1, 2, 1)', 'smd6mp(1, 0, 2, 1)', };

% tic;
blmapping_trueEvaldemo('smd7mp(1, 2, 1)',  22, 'use_seeding', false, 'seeding_strategy', 1);
% toc;

% % % problems = { 'smd8mp(1, 1, 1)'};
% % --------
seeds = 12: 21;
ns = length(seeds);
np = length(problems);
% ------

% create parameter for each method
paras  = struct([]);

% create parameter for the first method
for i = 1 : np
    for j = 1: ns
        strc_id = (i-1) * ns + j;
        paras(strc_id). problem_str = problems{i};
        paras(strc_id). seed = seeds(j);
        paras(strc_id). use_seeding = false;
        paras(strc_id).seeding_strategy = 0;
    end
end

for i = 1 : np
    for j = 1: ns
        strc_id = (i-1) * ns + j + ns * np * 1;
        paras(strc_id). problem_str = problems{i};
        paras(strc_id). seed = seeds(j); 
        paras(strc_id). use_seeding = true;
        paras(strc_id).seeding_strategy = 1;
    end
end


% create parameter for the second method
for i = 1 : np
    for j = 1: ns
        strc_id = (i-1) * ns + j + ns *np *2;
        paras(strc_id). problem_str = problems{i};
        paras(strc_id). seed = seeds(j);
        paras(strc_id). use_seeding = true;       % use local search
        paras(strc_id).seeding_strategy = 2;
    end
end


nrun = length(paras);
% 
tic;
parfor i = 1:nrun
    blmapping_trueEvaldemo(paras(i).problem_str, paras(i).seed,...
                             'use_seeding', paras(i).use_seeding, 'seeding_strategy',  paras(i).seeding_strategy);
                        
end
toc;




problems = {'smd5mp(1, 2, 1)' , 'smd7mp(1, 2, 1)',  'smd8mp(1, 2, 1)', ...
    'smd1mp(1, 2 , 1)' , 'smd2mp(1, 2, 1)',  'smd3mp(1, 2, 1)', ...
    'smd4mp(1, 2, 1)', 'smd6mp(1, 0, 2, 1)', };
% % --------

seeds = 12:21;
ns = length(seeds);
np = length(problems);
% ------

% create parameter for each method
paras  = struct([]);

% create parameter for the first method
for i = 1 : np
    for j = 1: ns
        strc_id = (i-1) * ns + j;
        paras(strc_id). problem_str = problems{i};
        paras(strc_id). seed = seeds(j);
        paras(strc_id). use_seeding = false;
        paras(strc_id).seeding_strategy = 0;
    end
end

for i = 1 : np
    for j = 1: ns
        strc_id = (i-1) * ns + j + ns * np * 1;
        paras(strc_id). problem_str = problems{i};
        paras(strc_id). seed = seeds(j); 
        paras(strc_id). use_seeding = true;
        paras(strc_id).seeding_strategy = 1;
    end
end


% create parameter for the second method
for i = 1 : np
    for j = 1: ns
        strc_id = (i-1) * ns + j + ns *np *2;
        paras(strc_id). problem_str = problems{i};
        paras(strc_id). seed = seeds(j);
        paras(strc_id). use_seeding = true;       % use local search
        paras(strc_id).seeding_strategy = 2;
    end
end


nrun = length(paras);
% 
tic;
parfor i = 1:nrun
    blmapping_trueEvaldemo(paras(i).problem_str, paras(i).seed,...
                             'use_seeding', paras(i).use_seeding, 'seeding_strategy',  paras(i).seeding_strategy);
                        
end
toc;


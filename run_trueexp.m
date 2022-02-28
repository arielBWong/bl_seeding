%%this script is to run ulego with multiple seeds
%
clear;
close all;

problem_folder = strcat(pwd, '/problems/TP3');
addpath(problem_folder);

problem_folder = strcat(pwd, '/problems/SMD');
addpath(problem_folder);

problem_folder = strcat(pwd, '/Surrogate/Methods/Surrogate');
addpath(problem_folder);

problem_folder = strcat(pwd, '/globalsolver');
addpath(problem_folder);

problem_folder = strcat(pwd, '/egosolver');
addpath(problem_folder);

problem_folder = strcat(pwd, '/ND_Sort');
addpath(problem_folder);

problem_folder = strcat(pwd,'/Utility');
addpath(problem_folder);

problem_folder = strcat(pwd,'/OODACE');
addpath(problem_folder);

problem_folder = strcat(pwd,'/transfer_utilities');
addpath(problem_folder);

%
% problems =   {'smd5(1, 2, 1)' , 'smd7(1, 2, 1)',  'smd8(1, 2, 1)', ...
%     'smd1(1, 2 , 1)' , 'smd2(1, 2, 1)',  'smd3(1, 2, 1)', ...
%     'smd4(1, 2, 1)', 'smd6(1, 0, 2, 1)', };

% tic;
% blmapping_egoEvaldemo('smd7mp(1, 1, 1)',  1, 'use_seeding', false, 'seeding_strategy', 3, 'threshold', 0.9);
% toc;
%



problem_sets = cell(1, 1);


problem_sets{1} = {'smd1mp(1, 1, 1)' , 'smd2mp(1, 1, 1)',  'smd3mp(1, 1, 1)', 'smd4mp(1, 1, 1)', ...
    'smd5mp(1, 1, 1)', 'smd6mp(1, 0, 1, 1)', 'smd7mp(1, 1, 1)',  'smd8mp(1, 1, 1)', ...
    'smd1(1, 1, 1)' , 'smd2(1, 1, 1)',  'smd3(1, 1, 1)', 'smd4(1, 1, 1)', 'smd5(1, 1, 1)',...
    'smd6(1, 0, 1, 1)', 'smd8(1, 1, 1)'};


problem_sets{1} = {'smd1(1, 2, 1)', 'smd2(1, 2, 1)',  'smd3(1, 2, 1)', ...
     'smd4(1, 2, 1)', 'smd5(1, 1, 1)', 'smd6(1, 0, 2, 1)', 'smd7(1, 2, 1)',  'smd8(1, 2, 1)',...
     'smd1mp(1, 1, 1)' , 'smd2mp(1, 1, 1)',  'smd3mp(1, 1, 1)', 'smd4mp(1, 1, 1)', ...
     'smd5mp(1, 1, 1)', 'smd6mp(1, 0, 1, 1)', 'smd7mp(1, 1, 1)',  'smd8mp(1, 1, 1)'};



thresholds = [-1.1, 0, 0.9];
thresholds = [-1.1];
tic;
for its = 1:length(thresholds)
    for ip = 1:length(problem_sets)
        problems = problem_sets{ip};
        seeds = 1:11;
        ns = length(seeds);
        np = length(problems);

        %reate parameter for each method
        paras  = struct([]);

        % method for cokrg
        for i = 1 : np
            for j = 1: ns
                strc_id = (i-1) * ns + j;
                paras(strc_id).problem_str = problems{i};
                paras(strc_id).seed = seeds(j);
                paras(strc_id).use_seeding = true;
                paras(strc_id).seeding_strategy = 2;
                paras(strc_id).threshold = thresholds(its);
            end
        end

        % method close optimal
        for i = 1 : np
            for j = 1: ns
                strc_id = (i-1) * ns + j + ns * np * 1;
                paras(strc_id). problem_str = problems{i};
                paras(strc_id). seed = seeds(j);
                paras(strc_id). use_seeding = true;
                paras(strc_id).seeding_strategy = 3;
                paras(strc_id).threshold = thresholds(its);
            end
        end


        nrun = length(paras);
        parfor i = 1:nrun
            blmapping_egoEvaldemo(paras(i).problem_str, paras(i).seed,...
                'use_seeding', paras(i).use_seeding, 'seeding_strategy',  paras(i).seeding_strategy,...
                'threshold', paras(i).threshold);

        end
    end
end

for ip = 1:length(problem_sets)
    problems = problem_sets{ip};
    seeds = 1:11;
    ns = length(seeds);
    np = length(problems);

    %reate parameter for each method
    paras  = struct([]);

    % method for cokrg
    for i = 1 : np
        for j = 1: ns
            strc_id = (i-1) * ns + j;
            paras(strc_id). problem_str = problems{i};
            paras(strc_id). seed = seeds(j);
            paras(strc_id). use_seeding = false;
            paras(strc_id).seeding_strategy = 5;
            paras(strc_id).threshold = 0;
        end
    end


    nrun = length(paras);
    parfor i = 1:nrun
        blmapping_egoEvaldemo(paras(i).problem_str, paras(i).seed,...
            'use_seeding', paras(i).use_seeding, 'seeding_strategy',  paras(i).seeding_strategy,...
            'threshold', paras(i).threshold);

    end
end
toc;




%%this script is to run ulego with multiple seeds
%
clear;
close all;

problem_folder = strcat(pwd,'/problems/TP3');
addpath(problem_folder);

problem_folder = strcat(pwd,'/problems/SMD');
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
%
% problems =   {'smd5(1, 2, 1)' , 'smd7(1, 2, 1)',  'smd8(1, 2, 1)', ...
%     'smd1(1, 2 , 1)' , 'smd2(1, 2, 1)',  'smd3(1, 2, 1)', ...
%     'smd4(1, 2, 1)', 'smd6(1, 0, 2, 1)', };

% tic;
% blmapping_trueEvaldemo('smd7(1, 1, 1)',  1, 'use_seeding', true, 'seeding_strategy', 2, 'threshold', 0);
% toc;


% problem_sets = cell(1, 1);


% problem_sets{1} =  {'smd5(1, 1, 1)' , 'smd7(1, 1, 1)',  'smd8(1, 1, 1)', ...
%     'smd1(1, 1, 1)' , 'smd2(1, 1, 1)',  'smd3(1, 1, 1)', ...
%     'smd4(1, 1, 1)', 'smd6(1, 0, 1, 1)', };

% 
% 
% thresholds = [0];
% tic;
% for its = 1:length(thresholds)
%     for ip = 1:length(problem_sets)
%         problems = problem_sets{ip};
%         seeds = 1:21;
%         ns = length(seeds);
%         np = length(problems);
% 
%         % create parameter for each method
%         paras  = struct([]);
% 
%         % method for cokrg
%         for i = 1 : np
%             for j = 1: ns
%                 strc_id = (i-1) * ns + j;
%                 paras(strc_id). problem_str = problems{i};
%                 paras(strc_id). seed = seeds(j);
%                 paras(strc_id). use_seeding = true;
%                 paras(strc_id).seeding_strategy = 1;
%                 paras(strc_id).threshold = thresholds(its);
%             end
%         end
% 
%         % method close optimal
%         for i = 1 : np
%             for j = 1: ns
%                 strc_id = (i-1) * ns + j + ns * np * 1;
%                 paras(strc_id). problem_str = problems{i};
%                 paras(strc_id). seed = seeds(j);
%                 paras(strc_id). use_seeding = false;
%                 paras(strc_id).seeding_strategy = 0;
%                 paras(strc_id).threshold = thresholds(its);
%             end
%         end
% 
% 
%         nrun = length(paras);
%         parfor i = 1:nrun
%             blmapping_trueEvaldemo(paras(i).problem_str, paras(i).seed,...
%                 'use_seeding', paras(i).use_seeding, 'seeding_strategy',  paras(i).seeding_strategy,...
%                 'threshold', paras(i).threshold);
% 
%         end
%     end
% end
% toc;


problem_sets = cell(1, 1);
problem_sets{1} =  {'smd5(1, 1, 1)' , 'smd7(1, 1, 1)',  'smd8(1, 1, 1)', ...
    'smd1(1, 1, 1)' , 'smd2(1, 1, 1)',  'smd3(1, 1, 1)', ...
    'smd4(1, 1, 1)', 'smd6(1, 0, 1, 1)', };

thresholds = [0, 0.1, 0.3, 0.5, 0.7, 0.9];
tic;
for its = 1:length(thresholds)
    for ip = 1:length(problem_sets)
        problems = problem_sets{ip};
        seeds = 1:21;
        ns = length(seeds);
        np = length(problems);

        % create parameter for each method
        paras  = struct([]);

        % method for cokrg
        for i = 1 : np
            for j = 1: ns
                strc_id = (i-1) * ns + j;
                paras(strc_id). problem_str = problems{i};
                paras(strc_id). seed = seeds(j);
                paras(strc_id). use_seeding = true;
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
            blmapping_trueEvaldemo(paras(i).problem_str, paras(i).seed,...
                'use_seeding', paras(i).use_seeding, 'seeding_strategy',  paras(i).seeding_strategy,...
                'threshold', paras(i).threshold);

        end
    end
end
toc;




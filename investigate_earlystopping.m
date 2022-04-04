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

seed = 1;
rng(seed, 'twister');

prob = smd1mp(1, 2, 1);
n = 20;
xu_pop  = repmat(prob.xu_bl, n, 1) + repmat(prob.xu_bu - prob.xu_bl, n, 1) .* lhsdesign(n, prob.n_uvar);


for i = 1:n
    xui = xu_pop(i, :);
    [match_xl, ~, ~, ~, ~] = llmatch_egoEvaluation(xui, prob, ...    
        'seeding_only', false,  'seeding_strategy', 7);

    xl_prime = prob.get_xlprime(xui);
    save_2xl = [xl_prime; match_xl];



    xl_prime = (xl_prime - prob.xl_bl) ./ (prob.xl_bu - prob.xl_bl);
    match_xl = (match_xl - prob.xl_bl) ./ (prob.xl_bu - prob.xl_bl);

    d = pdist2(xl_prime, match_xl);
    save_results(save_2xl, d, xui, prob, i);
   
end


function save_results(save_2xl, d, xui, prob, id)


name = strcat('investigate_earlystopping_nv_', num2str(prob.n_lvar));
resultfolder = fullfile(pwd, name );


n = exist(resultfolder);
if n ~= 7
    mkdir(resultfolder)
end

filename = strcat('xu_no_', num2str(id), '.csv');
savename = fullfile(resultfolder, filename);
csvwrite(savename, xui);

filename = strcat('save2xl_no_', num2str(id), '.csv');
savename = fullfile(resultfolder, filename);
csvwrite(savename, save_2xl);


filename = strcat('d_no_', num2str(id), '.csv');
savename = fullfile(resultfolder, filename);
csvwrite(savename, d);


end


% this m file is used to generate an example figure of cec2022
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

problem_folder = strcat(pwd,'/transfer_utilities');
addpath(problem_folder);


blmapping_trueEvaldemo('smd7mp(1, 1, 1)',  5, 'use_seeding', true, 'seeding_strategy', 3, 'threshold', 0.5);
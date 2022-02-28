% generate accuracy analysis
problem_folder = strcat(pwd,'/problems/TP3');
addpath(problem_folder);


thr = -110;
problems_2 = cell(1, 1);

% problems_2{1} =  { 'smd1(1, 2 , 1)' , 'smd2(1, 2, 1)',  'smd3(1, 2, 1)', ...
%      'smd4(1, 2, 1)','smd5(1, 2, 1)' , 'smd6(1, 0, 2, 1)', 'smd7(1, 2, 1)', 'smd8(1, 2, 1)' };

% problems_2{1} =  {'smd1(1, 1, 1)' , 'smd2(1, 1, 1)',  'smd3(1, 1, 1)', ...
%     'smd4(1, 1, 1)', 'smd5(1, 1, 1)' , 'smd6(1, 0, 1, 1)',  'smd7(1, 1, 1)',  'smd8(1, 1, 1)'};


% problems_2{1} = {'smd1(1, 2 , 1)' , 'smd2(1, 2, 1)',  'smd3(1, 2, 1)', ...
%                   'smd4(1, 2, 1)','smd5(1, 2, 1)' , 'smd6(1, 0, 2, 1)', 'smd7(1, 2, 1)', 'smd8(1, 2, 1)',...
%                   'smd1mp(1, 2, 1)' , 'smd2mp(1, 2, 1)',  'smd3mp(1, 2, 1)',  'smd4mp(1, 2, 1)', ....
%                  'smd5mp(1, 2, 1)' , 'smd6mp(1, 0, 2, 1)', 'smd7mp(1, 2, 1)',  'smd8mp(1, 2, 1)'};


% problems_2{1} = { 'smd1(1, 1, 1)' , 'smd2(1, 1, 1)',  'smd3(1, 1, 1)', ...
%      'smd4(1, 1, 1)', 'smd5(1, 1, 1)' , 'smd6(1, 0, 1, 1)',  'smd7(1, 1, 1)', ...
%      'smd8(1, 1, 1)',...
%     'smd1mp(1, 1, 1)' , 'smd2mp(1, 1, 1)',  'smd3mp(1, 1, 1)',  'smd4mp(1, 1, 1)', ....
%     'smd5mp(1, 1, 1)', 'smd6mp(1, 0, 1, 1)', 'smd7mp(1, 1, 1)',  'smd8mp(1, 1, 1)'};
%
%
problems_2{1} = {'smd1mp(1, 1, 1)' , 'smd2mp(1, 1, 1)',  'smd3mp(1, 1, 1)',  'smd4mp(1, 1, 1)', ....
    'smd5mp(1, 1, 1)', 'smd6mp(1, 0, 1, 1)', 'smd7mp(1, 1, 1)',  'smd8mp(1, 1, 1)'};

problems_2{1} = { 'smd1mp(1, 1, 1)' , 'smd2mp(1, 1, 1)',  'smd3mp(1, 1, 1)', 'smd4mp(1, 1, 1)', ...
    'smd5mp(1, 1, 1)', 'smd6mp(1, 0, 1, 1)', 'smd7mp(1, 1, 1)',  'smd8mp(1, 1, 1)', ...
    'smd1(1, 1, 1)' , 'smd2(1, 1, 1)',  'smd3(1, 1, 1)', 'smd4(1, 1, 1)', 'smd5(1, 1, 1)',...
    'smd6(1, 0, 1, 1)', 'smd7(1, 1, 1)', 'smd8(1, 1, 1)'};

% problems_2{1} = {'smd1mp(1, 2, 1)' , 'smd2mp(1, 2, 1)',  'smd3mp(1, 2, 1)',  'smd4mp(1, 2, 1)', ....
%               'smd5mp(1, 2, 1)' , 'smd6mp(1, 0, 2, 1)', 'smd7mp(1, 2, 1)',  'smd8mp(1, 2, 1)'};
% methods = {'_baseline_ea', '_seeding_strategy_1', '_seeding_strategy_2', '_seeding_strategy_3'};

for is = 1:length(problems_2)
    % methods = {'_baseline_ea', '_seeding_strategy_1', '_seeding_strategy_3'};
    methods = {'_baseline_ea', '_seeding_strategy_2', '_seeding_strategy_3'};
    problems = problems_2{is};
    prob_test = eval(problems{1});
    nv = prob_test.n_lvar;

    np = length(problems);
    seed = 11;
    mseed = 5;
    sigTestIndex = 3;  % refer to the newest algorithm which is 4 in this case


    foldername = strcat('resultfolder_trueEval', num2str(nv), '_thr_', num2str(thr));
    resultfolder = fullfile(pwd, foldername);

    outfoldername = strcat('processedresult_trueEval', num2str(nv),'_thr_', num2str(thr));
    outfoldername = fullfile(pwd, outfoldername);

    if ~exist(outfoldername, 'dir')
        mkdir(outfoldername);
    end

    % accuracy_extraction(problems, methods, resultfolder, np, seed, mseed, sigTestIndex, outfoldername);

    % methods = { '_seeding_strategy_2',  '_seeding_strategy_3'};
    % methods = {  '_seeding_strategy_3'};
    FE_analysis(problems, methods, resultfolder, np, seed, mseed, outfoldername) ;

    % methods = { '_seeding_strategy_2', '_seeding_strategy_3', };
    % switch_ratio(problems, methods, resultfolder, seed, mseed,outfoldername);

    % methods = {'_baseline_ea',  '_seeding_strategy_1'};
    % acccuracy_baseline(problems, methods, resultfolder, np, seed, mseed, 2, outfoldername);
end



function[] = acccuracy_baseline(problems, method, resultfolder, np, seed, mseed, sigTestIndex,outfoldername)
% sigTestIndex is to specify sigtest columne,


nm = length(method);
permethod_accuracy_up = {}; % upper save
permethod_accuracy_down = {}; % lower save

prob = eval(problems{1});
prefix = 'baseline';

for m = 1:nm
    accuracy_up  = zeros(np, seed); % for one problem [problem, seed]
    accuracy_low = zeros(np, seed);

    for p = 1:np
        prob = eval(problems{p});
        for s = 1: seed
            filename = strcat('final_accuracy_seed_', num2str(s), '.csv');
            foldername =  strcat(prob.name, method{m});
            savename = fullfile(resultfolder, foldername, filename);
            accuracy = csvread(savename);

            accuracy_up(p, s) = accuracy(1);
            accuracy_low(p, s) = accuracy(2);
        end
    end

    permethod_accuracy_up{m} = accuracy_up;
    permethod_accuracy_down{m} = accuracy_low;
end



filename = strcat(prefix, 'baselineMedian_accuracy_latex_nlvar_', num2str(prob.n_lvar),'.csv');
% filename = fullfile(outfoldername, filename);
fp = fopen(filename, 'w');
fprintf(fp, 'Problems,  baseline,  , neighbour optimal,  , sig\n');

for i = 1:np
    prob = eval(problems{i});
    fprintf(fp, '%s &,', prob.name);

    for j = 1:nm

        if j == 1
            [~, id] = sort(permethod_accuracy_up{j}(i, :));
            nx = id(mseed);
            ula1 =  permethod_accuracy_up{j}(i, nx);
            lla1 =  permethod_accuracy_down{j}(i, nx);

            % only suits 2 methods situation
            [~, id] = sort(permethod_accuracy_up{j+1}(i, :));
            nx = id(mseed);
            ula2 =  permethod_accuracy_up{j+1}(i, nx);
            lla2 =  permethod_accuracy_down{j+1}(i, nx);

            if ula1 < ula2
                fprintf(fp, '\\textbf{%0.4f} & , ', ula1);
            else
                fprintf(fp, '%0.4f & , ', ula1);
            end

            if lla1 < lla2
                fprintf(fp, '\\textbf{%0.4f} & , ', lla1);
            else
                fprintf(fp, '%0.4f & , ', lla1);
            end


        end

        if j == 2
            [~, id] = sort(permethod_accuracy_up{j}(i, :));
            nx = id(mseed);
            ula2 =  permethod_accuracy_up{j}(i, nx);
            lla2 =  permethod_accuracy_down{j}(i, nx);



            [~, id] = sort(permethod_accuracy_up{j-1}(i, :));
            nx = id(mseed);
            ula1 =  permethod_accuracy_up{j-1}(i, nx);
            lla1 =  permethod_accuracy_down{j-1}(i, nx);

            if ula2 < ula1
                fprintf(fp, '\\textbf{%0.4f} & , ', ula2);
            else
                fprintf(fp, '%0.4f & , ', ula2);
            end

            if lla2 < lla1
                fprintf(fp, '\\textbf{%0.4f} & , ', lla2);
            else
                fprintf(fp, '%0.4f & , ', lla2);
            end

        end


        if j == sigTestIndex  % accuracy smaller the better
            % 1 means  the former is smaller,  second is bigger/better,
            % for hv comparison larger the better
            [p1, h1, stats1] = ranksum( permethod_accuracy_up{j}(i, :), permethod_accuracy_up{j-1}(i, :),  'alpha', 0.05, 'tail', 'left');
            [p2, h2, stats2] = ranksum( permethod_accuracy_up{j-1}(i, :), permethod_accuracy_up{j}(i, :),  'alpha', 0.05, 'tail', 'left');
            if h1 == 1 && h2 == 0
                fprintf(fp, '$\\uparrow$, ');
            elseif h2==1 && h1 == 0
                fprintf(fp, '$\\downarrow$,');
            else
                fprintf(fp, '$\\approx$, ');
            end
        end


    end
    fprintf(fp, '\\\\\n');

end
fclose(fp);


% normalized accuracy
%
norm_up = [];
norm_down = [];
for i = 1:nm
    norm_up = [norm_up; permethod_accuracy_up{i}]; %[problem, seed]
    norm_down = [norm_down; permethod_accuracy_down{i}];
end
norm_up = (norm_up - min(norm_up(:)))/(max(norm_up(:)) - min(norm_up(:)));
norm_down = (norm_down - min(norm_down(:)))/(max(norm_down(:)) - min(norm_down(:)));
norm_accuracy = norm_up + norm_down;

permethod_normAccuracy = {};

for i = 1:nm
    permethod_normAccuracy{i} = norm_accuracy((i-1)*np + 1: i*np, :);
end

filename = strcat(prefix, 'median_accuracy_Normlatex_', num2str(prob.n_lvar),'.csv');
filename = fullfile(outfoldername, filename);
fp = fopen(filename, 'w');
fprintf(fp, 'problems,  baseline,  neighbour optimal, decision and cokring, decision and neighbour optimal, \n');
for i = 1:np
    prob = eval(problems{i});
    fprintf(fp, '%s &,', prob.name);

    for j = 1:nm

        [~, id] = sort(permethod_normAccuracy{j}(i, :));
        nx = id(mseed);
        fprintf(fp, '%0.4f & , ', permethod_normAccuracy{j}(i, nx));

        if j == sigTestIndex %
            %  1 means  the former is smaller,  second is bigger/better,
            % for accuracy smaller better
            [p1,h1,stats1] = ranksum( permethod_normAccuracy{j}(i, :), permethod_normAccuracy{j-1}(i, :),  'alpha', 0.05, 'tail', 'left');
            [p2,h2,stats2] = ranksum( permethod_normAccuracy{j-1}(i, :), permethod_normAccuracy{j}(i, :),  'alpha', 0.05, 'tail', 'left');
            if h1 == 1 && h2 == 0
                fprintf(fp, '%s $\\uparrow$ %s, ', 'DMNE', 'DMCK');
            elseif h2==1 && h1 == 0
                fprintf(fp, '%s $\\downarrow$ %s, ', 'DMNE', 'DMCK');
            else
                fprintf(fp, '%s $\\approx$ %s, ', 'DMNE', 'DMCK');
            end
        end
    end
    fprintf(fp, '\\\\\n');

end
fclose(fp);

end

function[] = switch_ratio(problems, methods, resultfolder, ns, median_ns, outfoldername)
% This function calculates the ratio of switch to global search after
% the first generation
% only work for strategy 2 and 3
%----------

prob = eval(problems{1});
diff = 'mp';
if contains(prob.name, diff)
    prefix = 'smd_mp_';
else
    prefix = 'smd_';
end


np = length(problems);
nm = length(methods);
% cell(num_problems):  [num_seeds, num_methods]
raw_switch = cell(1, np);
% cell(num_problems):[1, num_methods];
std_switch = cell(1, np);
% cell(num_problems):  [1, num_methods];
mean_switch = cell(1, np);


for ip = 1:np

    prob = eval(problems{ip});
    raw_switch_single = zeros(ns, nm);

    for is = 1:ns
        for im = 1:nm
            foldername= fullfile(resultfolder, strcat(prob.name, methods{im}));
            filename = fullfile(foldername, strcat('lowerlevelswitch_seed_', num2str(is), '.mat'));
            load(filename);
            % for variable lower_decisionSwitch 1 means that
            total_global = sum(lower_decisionSwitch(:));
            raw_switch_single(is, im) = total_global;
            clear lower_decisionSwitch;
        end
    end

    raw_switch{ip} = raw_switch_single;
    std_switch{ip} = std(raw_switch_single, 1);
    mean_switch{ip} = mean(raw_switch_single, 1);
end

% write into files
filename = strcat(prefix, 'globalSwitch_mean_nvar_', num2str(prob.n_lvar),'.csv');
filename = fullfile(outfoldername, filename);
fp = fopen(filename, 'w');
fprintf(fp, 'problem/methods, ');
for im = 1:nm
    fprintf(fp,  '%s, ', methods{im});
end
fprintf(fp, '\n');

for ip = 1:np
    prob = eval(problems{ip});
    fprintf(fp, '%s,',  prob.name);
    for im = 1:nm
        fprintf(fp, '%f, ',  floor(mean_switch{ip}(im)));
    end
    fprintf(fp, '\n');
end
fclose(fp);

end

function [] = lowerSuccessRate(problems, method, resultfolder, np, seed, mseed)

localsearch_success1 = {};
localsearch_success2 = {};
nm = length(method);


for m = 1:nm

    lowerSC1 = zeros(np, seed); % number of local search invoked
    lowerSC2 = zeros(np, seed); % ratio of successful local search

    for p = 1:np
        prob = eval(problems{p});
        for s = 1: seed

            filename = strcat('lower_success_seed_', num2str(s), '.csv');
            foldername =  strcat(prob.name, method{m});
            savename = fullfile(resultfolder, foldername, filename);
            lower_success = csvread(savename);

            lowerSC1(p, s) = size(lower_success, 1);
            lowerSC2(p, s) = sum(lower_success)/size(lower_success, 1);

        end
    end

    localsearch_success1{m} = lowerSC1;
    localsearch_success2{m} = lowerSC2;
end


% analysis median number of local search, median value of local search
% success rate
filename = strcat('success_rate_check_nlvar_', num2str(prob.n_lvar),'.csv');

fp = fopen(filename, 'w');
fprintf(fp, 'problems,  transfer to local search number, rate, transfer to local search+ number, rate, compare,  decision making number, rate, compare\n');
for i = 1:np
    prob = eval(problems{i});
    fprintf(fp, '%s &,', prob.name);

    for j = 1:nm
        [~, id] = sort(localsearch_success1{j}(i, :)); %number
        nx = id(mseed);
        fprintf(fp, '%d & ,', localsearch_success1{j}(i, nx));

        [~, id] = sort(localsearch_success2{j}(i, :)); %rate
        nx = id(mseed);
        fprintf(fp, '%0.4f & ,', localsearch_success2{j}(i, nx));

        if j == 3
            % conduct significant test cu
            % 1 means  the former is smaller,  second is bigger/better,
            % for hv comparison larger the better
            [p1,h1,stats1] = ranksum( localsearch_success2{j}(i, :), localsearch_success2{j-1}(i, :),  'alpha', 0.05, 'tail', 'left');
            [p2,h2,stats2] = ranksum( localsearch_success2{j-1}(i, :), localsearch_success2{j}(i, :),  'alpha', 0.05, 'tail', 'left');
            if h1 == 1 && h2 == 0
                fprintf(fp, '$\\downarrow_{4}$,');
            elseif h2==1 && h1 == 0
                fprintf(fp, '$\\uparrow_{4}$,');
            else
                fprintf(fp, '$\\approx_{4}$,');
            end


        end

        if j == 2
            % conduct significant test cu
            % 1 means  the former is smaller,  second is bigger/better,
            % for hv comparison larger the better
            [p1,h1,stats1] = ranksum( localsearch_success2{j}(i, :), localsearch_success2{j-1}(i, :),  'alpha', 0.05, 'tail', 'left');
            [p2,h2,stats2] = ranksum( localsearch_success2{j-1}(i, :), localsearch_success2{j}(i, :),  'alpha', 0.05, 'tail', 'left');
            if h1 == 1 && h2 == 0
                fprintf(fp, '$\\downarrow_{3}$ &,');
            elseif h2==1 && h1 == 0
                fprintf(fp, '$\\uparrow_{3}$ & ,');
            else
                fprintf(fp, '$\\approx_{3}$ & ,');
            end


        end

    end
    fprintf(fp, '\\\\\n');
end
fclose(fp);

end

function [] = lowerSuccessRateExtension(problems, method, resultfolder, np, seed, mseed)
localsearch_success1 = {};
localsearch_success2 = {};
nm = length(method);


for m = 1:nm

    lowerSC1 = zeros(np, seed); % number of local search invoked
    lowerSC2 = zeros(np, seed); % ratio of successful local search

    for p = 1:np
        prob = eval(problems{p});
        for s = 1: seed

            filename = strcat('lower_success_seed_', num2str(s), '.csv');
            foldername =  strcat(prob.name, method{m});
            savename = fullfile(resultfolder, foldername, filename);
            lower_success = csvread(savename);

            lowerSC1(p, s) = size(lower_success, 1);
            lowerSC2(p, s) = sum(lower_success)/size(lower_success, 1);

        end
    end

    localsearch_success1{m} = lowerSC1;
    localsearch_success2{m} = lowerSC2;
end


% analysis median number of local search, median value of local search
% success rate
filename = strcat('success_rate_check_nlvar_', num2str(prob.n_lvar),'.csv');

fp = fopen(filename, 'w');
fprintf(fp, 'problems,  transfer to local search number, rate, transfer to local search+ number, rate, compare,  decision making number, rate, compare, correction making number, rate, compare\n');
for i = 1:np
    prob = eval(problems{i});
    fprintf(fp, '%s &,', prob.name);

    for j = 1:nm
        [~, id] = sort(localsearch_success1{j}(i, :)); %number
        nx = id(mseed);
        fprintf(fp, '%d & ,', localsearch_success1{j}(i, nx));

        [~, id] = sort(localsearch_success2{j}(i, :)); %rate
        nx = id(mseed);
        fprintf(fp, '%0.4f & ,', localsearch_success2{j}(i, nx));

        if j == 3
            % conduct significant test cu
            % 1 means  the former is smaller,  second is bigger/better,
            % for hv comparison larger the better
            [p1,h1,stats1] = ranksum( localsearch_success2{j}(i, :), localsearch_success2{j-1}(i, :),  'alpha', 0.05, 'tail', 'left');
            [p2,h2,stats2] = ranksum( localsearch_success2{j-1}(i, :), localsearch_success2{j}(i, :),  'alpha', 0.05, 'tail', 'left');
            if h1 == 1 && h2 == 0
                fprintf(fp, '$\\downarrow_{4}$ &,');
            elseif h2==1 && h1 == 0
                fprintf(fp, '$\\uparrow_{4}$ &,');
            else
                fprintf(fp, '$\\approx_{4}$&,');
            end


        end

        if j == 2
            % conduct significant test cu
            % 1 means  the former is smaller,  second is bigger/better,
            % for hv comparison larger the better
            [p1,h1,stats1] = ranksum( localsearch_success2{j}(i, :), localsearch_success2{j-1}(i, :),  'alpha', 0.05, 'tail', 'left');
            [p2,h2,stats2] = ranksum( localsearch_success2{j-1}(i, :), localsearch_success2{j}(i, :),  'alpha', 0.05, 'tail', 'left');
            if h1 == 1 && h2 == 0
                fprintf(fp, '$\\downarrow_{3}$ &,');
            elseif h2==1 && h1 == 0
                fprintf(fp, '$\\uparrow_{3}$ & ,');
            else
                fprintf(fp, '$\\approx_{3}$ & ,');
            end


        end

        if j == 4
            % conduct significant test cu
            % 1 means  the former is smaller,  second is bigger/better,
            % for hv comparison larger the better
            [p1,h1,stats1] = ranksum( localsearch_success2{j}(i, :), localsearch_success2{j-1}(i, :),  'alpha', 0.05, 'tail', 'left');
            [p2,h2,stats2] = ranksum( localsearch_success2{j-1}(i, :), localsearch_success2{j}(i, :),  'alpha', 0.05, 'tail', 'left');
            if h1 == 1 && h2 == 0
                fprintf(fp, '$\\downarrow_{5}$,');
            elseif h2==1 && h1 == 0
                fprintf(fp, '$\\uparrow_{5}$,');
            else
                fprintf(fp, '$\\approx_{5}$,');
            end


        end

    end
    fprintf(fp, '\\\\\n');
end
fclose(fp);
end

function [] = FE_analysis(problems, method, resultfolder, np, seed, mseed, outfoldername)

prob = eval(problems{1});
% diff = 'mp';
% if contains(prob.name, diff)
%     prefix = 'smd_mp_';
% else
%     prefix = 'smd_';
% end

nm = length(method);
FE_results2 = {};
for m = 1:nm
    lowerFE  = zeros(np, seed);
    % lowerFE_extra = zeros(np, seed);

    for p = 1:np
        prob = eval(problems{p});
        for s = 1: seed
            filename = strcat('lowerlevelcount_seed_', num2str(s), '.csv');
            foldername =  strcat(prob.name, method{m});
            savename = fullfile(resultfolder, foldername, filename);
            two_counts = csvread(savename);
            lowerFE(p, s) = sum(two_counts);
            % lowerFE_extra(p, s) = lowerFE(2);
        end
    end
    FE_results2{m} = lowerFE;
end

filename = strcat('median_FE_nlvar_', num2str(prob.n_lvar),'.csv');
filename = fullfile(outfoldername, filename);
fp = fopen(filename, 'w');
% fprintf(fp, 'problems,');
% for im = 1:nm
%     fprintf(fp, '%s,', method{im});
% end
% fprintf(fp, '\n');


fprintf(fp, 'problems, global krg,  Correlation Cokrg, Correlation Neighbor \n');

for i = 1:np
    prob = eval(problems{i});
    fprintf(fp, '%s &,', prob.name);

    [~, id] = sort(FE_results2{1}(i, :));
    nx = id(mseed);
    fe1 =  FE_results2{1}(i, nx);

    [~, id] = sort(FE_results2{2}(i, :));
    nx = id(mseed);
    fe2 =  FE_results2{2}(i, nx);

    [~, id] = sort(FE_results2{3}(i, :));
    nx = id(mseed);
    fe3 =  FE_results2{3}(i, nx);

    fe = [fe1, fe2, fe3];
    [fes, idx] = sort(fe);


    for j = 1:nm

        if idx(1) == j
            fprintf(fp, '\\textbf{%d}  & , ', fe(j));
        else
            fprintf(fp, '%d &, ', fe(j));
        end

    end
    fprintf(fp, '\\\\\n');
    % fprintf(fp, '\n');

end
fclose(fp);
end

function[] = accuracy_extraction(problems, method, resultfolder, np, seed, mseed, sigTestIndex,outfoldername)
% sigTestIndex is to specify sigtest columne,


nm = length(method);
permethod_accuracy_up = {}; % upper save
permethod_accuracy_down = {}; % lower save

prob = eval(problems{1});
diff = 'mp';
prefix = 'mix';

% if contains(prob.name, diff)
%     prefix = 'smd_mp_';
% else
%     prefix = 'smd_';
% end

for m = 1:nm
    accuracy_up = zeros(np, seed); % for one problem [problem, seed]
    accuracy_low = zeros(np, seed);

    for p = 1:np
        prob = eval(problems{p});
        for s = 1: seed
            filename = strcat('final_accuracy_seed_', num2str(s), '.csv');
            foldername = strcat(prob.name, method{m});
            savename = fullfile(resultfolder, foldername, filename);
            accuracy = csvread(savename);

            accuracy_up(p, s) = accuracy(1);
            accuracy_low(p, s) = accuracy(2);
        end
    end

    permethod_accuracy_up{m} = accuracy_up;
    permethod_accuracy_down{m} = accuracy_low;
end



filename = strcat(prefix, 'median_accuracy_latex_nlvar_', num2str(prob.n_lvar),'.csv');
filename = fullfile(outfoldername, filename);
fp = fopen(filename, 'w');
fprintf(fp, 'problems, global krg,  ,Correlation Cokrg, ,  sig(2-1),Correlation Neighbor, , sig(3-2) \n');

for i = 1:np
    prob = eval(problems{i});
    fprintf(fp, ' %s &,', prob.name);

    [~, id] = sort(permethod_accuracy_up{1}(i, :));
    nx = id(mseed);
    ula1 =  permethod_accuracy_up{1}(i, nx);
    lla1 =  permethod_accuracy_down{1}(i, nx);

    % only suits 2 methods situation
    [~, id] = sort(permethod_accuracy_up{2}(i, :));
    nx = id(mseed);
    ula2 =  permethod_accuracy_up{2}(i, nx);
    lla2 =  permethod_accuracy_down{2}(i, nx);

    % only suits 3  methods situation
    [~, id] = sort(permethod_accuracy_up{3}(i, :));
    nx = id(mseed);
    ula3 =  permethod_accuracy_up{3}(i, nx);
    lla3 =  permethod_accuracy_down{3}(i, nx);

    ula = [ula1, ula2, ula3];
    lla = [lla1, lla2, lla3];

    [~, idxul] = sort(ula);
    [~, idxll] = sort(lla);

    for j = 1:nm
        if j == 1
            if idxul(1) ==j
                fprintf(fp, '\\textbf{%0.4f} & , ', ula1);
            else
                fprintf(fp, '%0.4f & , ', ula1);
            end

            if idxll(1) ==j
                fprintf(fp, '\\textbf{%0.4f} & , ', lla1);
            else
                fprintf(fp, '%0.4f & , ', lla1);
            end
        end

        if j == 2
            if idxul(1) ==j
                fprintf(fp, '\\textbf{%0.4f} & , ', ula2);
            else
                fprintf(fp, '%0.4f & , ', ula2);
            end

            if idxll(1) ==j
                fprintf(fp, '\\textbf{%0.4f} & , ', lla2);
            else
                fprintf(fp, '%0.4f & , ', lla2);
            end

            % compare the second strategy to the first strategy
            [p1, h1, stats1] = ranksum( permethod_accuracy_up{j}(i, :), permethod_accuracy_up{j-1}(i, :),  'alpha', 0.05, 'tail', 'left');
            [p2, h2, stats2] = ranksum( permethod_accuracy_up{j-1}(i, :), permethod_accuracy_up{j}(i, :),  'alpha', 0.05, 'tail', 'left');
            if h1 == 1 && h2 == 0
                fprintf(fp, '$\\uparrow$ &, ');
            elseif h2==1 && h1 == 0
                fprintf(fp, '$\\downarrow$ &, ');
            else
                fprintf(fp, '$\\approx$ &, ');
            end

        end


        if j == 3
            if idxul(1) ==j
                fprintf(fp, '\\textbf{%0.4f} & , ', ula3);
            else
                fprintf(fp, '%0.4f & , ', ula3);
            end

            if idxll(1) ==j
                fprintf(fp, '\\textbf{%0.4f} & , ', lla3);
            else
                fprintf(fp, '%0.4f & , ', lla3);
            end
        end


        if j == sigTestIndex  % accuracy smaller the better
            % 1 means  the former is smaller,  second is bigger/better,
            % for hv comparison larger the better
            [p1, h1, stats1] = ranksum( permethod_accuracy_up{j}(i, :), permethod_accuracy_up{j-1}(i, :),  'alpha', 0.05, 'tail', 'left');
            [p2, h2, stats2] = ranksum( permethod_accuracy_up{j-1}(i, :), permethod_accuracy_up{j}(i, :),  'alpha', 0.05, 'tail', 'left');
            if h1 == 1 && h2 == 0
                fprintf(fp, '$\\uparrow$, ');
            elseif h2==1 && h1 == 0
                fprintf(fp, '$\\downarrow$, ');
            else
                fprintf(fp, '$\\approx$, ');
            end
        end



    end
    fprintf(fp, '\\\\\n');

end
fclose(fp);


% normalized accuracy
%
norm_up = [];
norm_down = [];
for i = 1:nm
    norm_up = [norm_up; permethod_accuracy_up{i}]; %[problem, seed]
    norm_down = [norm_down; permethod_accuracy_down{i}];
end
norm_up = (norm_up - min(norm_up(:)))/(max(norm_up(:)) - min(norm_up(:)));
norm_down = (norm_down - min(norm_down(:)))/(max(norm_down(:)) - min(norm_down(:)));
norm_accuracy = norm_up + norm_down;

permethod_normAccuracy = {};

for i = 1:nm
    permethod_normAccuracy{i} = norm_accuracy((i-1)*np + 1: i*np, :);
end

filename = strcat(prefix, 'median_accuracy_Normlatex_', num2str(prob.n_lvar),'.csv');
filename = fullfile(outfoldername, filename);
fp = fopen(filename, 'w');
fprintf(fp, 'problems,  baseline,  neighbour optimal, decision and cokring, decision and neighbour optimal, \n');
for i = 1:np
    prob = eval(problems{i});
    fprintf(fp, '%s &,', prob.name);

    for j = 1:nm

        [~, id] = sort(permethod_normAccuracy{j}(i, :));
        nx = id(mseed);
        fprintf(fp, '%0.4f & , ', permethod_normAccuracy{j}(i, nx));

        if j == sigTestIndex %
            %  1 means  the former is smaller,  second is bigger/better,
            % for accuracy smaller better
            [p1,h1,stats1] = ranksum( permethod_normAccuracy{j}(i, :), permethod_normAccuracy{j-1}(i, :),  'alpha', 0.05, 'tail', 'left');
            [p2,h2,stats2] = ranksum( permethod_normAccuracy{j-1}(i, :), permethod_normAccuracy{j}(i, :),  'alpha', 0.05, 'tail', 'left');
            if h1 == 1 && h2 == 0
                fprintf(fp, '%s $\\uparrow$ %s, ', 'DMNE', 'DMCK');
            elseif h2==1 && h1 == 0
                fprintf(fp, '%s $\\downarrow$ %s, ', 'DMNE', 'DMCK');
            else
                fprintf(fp, '%s $\\approx$ %s, ', 'DMNE', 'DMCK');
            end
        end
    end
    fprintf(fp, '\\\\\n');

end
fclose(fp);



end

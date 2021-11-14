% generate accuracy analysis
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

% 
% problems = { 'smd1mp(1, 2, 1)' , 'smd2mp(1, 2, 1)',  'smd3mp(1, 2, 1)',  'smd4mp(1, 2, 1)', ....
%     'smd5mp(1, 2, 1)' , 'smd6mp(1, 0, 2, 1)', 'smd7mp(1, 2, 1)',  'smd8mp(1, 2, 1)'};


problems = { 'smd1mp(1, 1, 1)' , 'smd2mp(1, 1, 1)',  'smd3mp(1, 1, 1)',  'smd4mp(1, 1, 1)', ....
     'smd5mp(1, 1, 1)' , 'smd6mp(1, 0, 1, 1)', 'smd7mp(1, 1, 1)',  'smd8mp(1, 1, 1)'};

% problems = { 'smd1mp(1, 1, 1)'};
method = {'_seeding_1_rand_0', '_seeding_1_rand_1', '_seeding_3_rand_0',   '_seeding_2_rand_1'};
prob_test = eval(problems{1});
nv = prob_test.n_lvar;

foldername = strcat('resultfolder_trueEval', num2str(nv));
resultfolder = fullfile(pwd, foldername);


np = length(problems);
seed = 29;
mseed = 15;
sigTestIndex = 4;  % refer to the newest algorithm which is 4 in this case 

internal_check(problems, seed, resultfolder);

% lowerSuccessRate2(problems, method, resultfolder, np, seed, mseed)
% check_2successrate(problems, seed, resultfolder);

% accuracy_extraction(problems, method, resultfolder, np, seed, mseed, sigTestIndex);
% accuracy_extractionExtension(problems, method, resultfolder, np, seed, mseed, sigTestIndex);
% FE_analysis(problems, method, resultfolder, np, seed, mseed, sigTestIndex) 


% success rate
% method = {'_seeding_1_rand_0', '_seeding_1_rand_1', '_seeding_3_rand_0', '_seeding_2_rand_0', '_seeding_2_rand_1'};
% lowerSuccessRate(problems, method, resultfolder, np, seed, mseed);
% lowerSuccessRateExtension(problems, method, resultfolder, np, seed, mseed);


function [] = internal_check(problems, seed, resultfolder)
method =  '_seeding_2_rand_1';
prob = eval(problems{1});
filename = strcat('success_rate_checkkrg2Rand4opt_nlvar_', num2str(prob.n_lvar),'.csv');
np = length(problems);
fp = fopen(filename, 'w');
fprintf(fp, 'problems, Trans success & random success (11), Trans sucess & random fail (10) ,  Trans fail * random sucess (01), Trans fail and random fail (00)\n');
for p = 1:np
    option11 = zeros(seed, 1);
    option10 = zeros(seed, 1);
    option01 = zeros(seed, 1);
    option11 = zeros(seed, 1);
    
    prob = eval(problems{p});
    fprintf(fp, '%s,', prob.name);
    
    for s = 1: seed
        filename = strcat('lower_success2_seed_', num2str(s), '.csv');
        foldername =  strcat(prob.name, method);
        savename = fullfile(resultfolder, foldername, filename);
        lower_success2 = csvread(savename);
        
        %
        option11(s, 1) = sum(lower_success2(:,1) .* lower_success2(:, 2));
        
         krgsuccess_idx = lower_success2(:, 1) == 1;
        randomfail_idx = lower_success2(:, 2) == 0;
       
        
        option10_idx = randomfail_idx .* krgsuccess_idx;
        option10(s, 1) = sum(option10_idx);
        
        krgfail_idx = lower_success2(:, 1) == 0;
        randomsucess_idx = lower_success2(:, 2) == 1;
        
        option01_idx = krgfail_idx .* randomsucess_idx;
        option01(s, 1) = sum(option01_idx);
        
        option00_idx = krgfail_idx .* randomfail_idx;
        option00(s, 1) = sum(option00_idx);
    end
    
    fprintf(fp, '%d, ', median(option11));
    fprintf(fp, '%d, ', median(option10));
    fprintf(fp, '%d, ', median(option01));
    fprintf(fp, '%d\n ', median(option00));
end
fclose(fp);
end


function [] = check_2successrate(problems, seed, resultfolder)
method =  '_seeding_2_rand_1';
prob = eval(problems{1});
filename = strcat('success_rate_checkKrg2Rand_nlvar_', num2str(prob.n_lvar),'.csv');
np = length(problems);
fp = fopen(filename, 'w');
fprintf(fp, 'problems, cokrg success, random success, significance transfer to random\n');
for p = 1:np
    transfer_success = zeros(1, seed);
    random_success = zeros(1, seed);
    prob = eval(problems{p});
    fprintf(fp, '%s,', prob.name);
    
    for s = 1: seed
        filename = strcat('lower_success2_seed_', num2str(s), '.csv');
        foldername =  strcat(prob.name, method);
        savename = fullfile(resultfolder, foldername, filename);
        lower_success2 = csvread(savename);
        
        %
        transfer_success(s) = sum(lower_success2(:, 1));
        random_success(s) = sum(lower_success2(:, 2));
    end
    
    fprintf(fp, '%d, ', median(transfer_success));
    fprintf(fp, '%d, ', median(random_success));
    
    % conduct significant test cu
    % 1 means  the former is smaller,  second is bigger/better,
    % for success comparison larger the better
    [p1,h1,stats1] = ranksum( transfer_success,  random_success,  'alpha', 0.05, 'tail', 'left');
    [p2,h2,stats2] = ranksum( random_success,  transfer_success,  'alpha', 0.05, 'tail', 'left');
    if h1 == 1 && h2 == 0
        fprintf(fp, '$\\downarrow$,');
    elseif h2==1 && h1 == 0
        fprintf(fp, '$\\uparrow$,');
    else
        fprintf(fp, '$\\approx$,');
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
fprintf(fp, 'problems,  transfer local , rate, transfer local + , rate,  random local, rate\n');
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
    end
     fprintf(fp, '\\\\\n');
end
fclose(fp);
end


function [] =  lowerSuccessRate2(problems, method, resultfolder, np, seed, mseed)
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
            lowerSC2(p, s) = sum(lower_success);
            
        end
    end    
    
    localsearch_success1{m} = lowerSC1;
    localsearch_success2{m} = lowerSC2;
end


% analysis median number of local search, median value of local search
% success rate
filename = strcat('success_rate_check24_nlvar_', num2str(prob.n_lvar),'.csv');

fp = fopen(filename, 'w');
fprintf(fp, 'problems,  trans local,  cokrg local, significance\n');
for i = 1:np
    prob = eval(problems{i});
    fprintf(fp, '%s &,', prob.name);
    
    j = 2;
    k = 4;
    
    [~, id] = sort(localsearch_success2{j}(i, :)); %number 
    nx = id(mseed);
    fprintf(fp, '%d,', localsearch_success2{j}(i, nx));
            
    [~, id] = sort(localsearch_success2{k}(i, :)); %rate 
    nx = id(mseed);
    fprintf(fp, '%d & ,', localsearch_success2{k}(i, nx));
 
        
     % conduct significant test cu
    % 1 means  the former is smaller,  second is bigger,
    % for success comparison larger the better
    [p1,h1,stats1] = ranksum(localsearch_success2{j}(i, :), localsearch_success2{k}(i, :),  'alpha', 0.05, 'tail', 'left');
    [p2,h2,stats2] = ranksum(localsearch_success2{k}(i, :), localsearch_success2{j}(i, :), 'alpha', 0.05, 'tail', 'left');
    
    if h1 == 1 && h2 == 0
        fprintf(fp, '$\\downarrow$,');
    elseif h2 == 1 && h1 == 0
        fprintf(fp, '$\\up$,');
    else
        fprintf(fp, '$\\approx$,');
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
filename = strcat('success_rate_check2_nlvar_', num2str(prob.n_lvar),'.csv');

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

function [] = FE_analysis(problems, method, resultfolder, np, seed, mseed, sigTestIndex) 

nm = length(method);

FE_results2 = {};
for m = 1:nm
    lowerFE  = zeros(np, seed);
   
    for p = 1:np
        prob = eval(problems{p});
        for s = 1: seed
            filename = strcat('lowerlevelcount_seed_', num2str(s), '.csv');
            foldername =  strcat(prob.name, method{m});
            savename = fullfile(resultfolder, foldername, filename);
            lowerFE(p, s) = csvread(savename);  
        end
    end
     FE_results2{m} = lowerFE;
end

filename = strcat('median_FE_nlvar_', num2str(prob.n_lvar),'.csv');
fp = fopen(filename, 'w');
fprintf(fp, 'problems, baseline,  transferred solution to local search, transferred solution to local search+ , selective transfer \n');
for i = 1:np
    prob = eval(problems{i});
    fprintf(fp, '%s &,', prob.name);
    
    for j = 1:nm
        
        [~, id] = sort(FE_results2{j}(i, :));
        nx = id(mseed);
        fprintf(fp, '%d &, ', FE_results2{j}(i, nx));
       
    end
    fprintf(fp, '\\\\\n');

end
fclose(fp);

end

function[] = accuracy_extraction(problems, method, resultfolder, np, seed, mseed, sigTestIndex)
% sigTestIndex is to specify sigtest columne,


nm = length(method);
method_results1 = {}; % upper save
method_results2 = {}; % lower save


for m = 1:nm
    accuracy_up  = zeros(np, seed);
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
    
    method_results1{m} = accuracy_up;
    method_results2{m} = accuracy_low;
    
end

filename = strcat('median_accuracy_latex_nlvar_', num2str(prob.n_lvar),'.csv');
fp = fopen(filename, 'w');
fprintf(fp, 'problems,  transfer to local search(Upper), (Lower), transfer to local search+ (Upper), (Lower), random to local (Upper), (Lower), Cokrg to local (Upper), (Lower)  \n');

for i = 1:np
    prob = eval(problems{i});
    fprintf(fp, '%s &,', prob.name);
    
    for j = 1:nm
        
        [~, id] = sort(method_results1{j}(i, :));
        nx = id(mseed);
        fprintf(fp, '%0.4f & , ',  method_results1{j}(i, nx));
        fprintf(fp, '%0.4f & , ',  method_results2{j}(i, nx));
 
    end
    fprintf(fp, '\\\\\n');

end
fclose(fp);

filename = strcat('median_accuracy_sigTest_random2Cokrg_nlvar_', num2str(prob.n_lvar),'.csv');
fp = fopen(filename, 'w');
fprintf(fp, 'problems,  random to local search(Upper), (Lower), Cokrg to local search (Upper), (Lower), Significance test \n');
% normalized accuracy
for i = 1:np
    prob = eval(problems{i});
    fprintf(fp, '%s,', prob.name);
    
    j = 1;
    [~, id] = sort(method_results1{j}(i, :));
    nx = id(mseed);
    fprintf(fp, '%0.4f & , ',  method_results1{j}(i, nx));
    fprintf(fp, '%0.4f & , ',  method_results2{j}(i, nx));
    
    j = 2;
    [~, id] = sort(method_results1{j}(i, :));
    nx = id(mseed);
    fprintf(fp, '%0.4f & , ',  method_results1{j}(i, nx));
    fprintf(fp, '%0.4f & , ',  method_results2{j}(i, nx));
    
    % conduct significant test cu
    % 1 means  the former is smaller,  second is bigger,
    % for accuracy comparison smaller the better
    [p1,h1,stats1] = ranksum( method_results1{3}(i, :), method_results1{4}(i, :),  'alpha', 0.05, 'tail', 'left');
    [p2,h2,stats2] = ranksum(method_results1{4}(i, :), method_results1{3}(i, :),  'alpha', 0.05, 'tail', 'left');
    if h1 == 1 && h2 == 0
        fprintf(fp, '$\\uparrow$,');
    elseif h2 == 1 && h1 == 0
        fprintf(fp, '$\\downarrow$,');
    else
        fprintf(fp, '$\\approx$,');
    end
    
    fprintf(fp, '\n');
end
fclose(fp);
%-----------------

% filename = strcat('median_accuracy_sigTest_Transfer2random_nlvar_', num2str(prob.n_lvar),'.csv');
% fp = fopen(filename, 'w');
% fprintf(fp, 'problems,  transfer to local search(Upper), (Lower), Random (Upper), (Lower), Significance test \n');
% % normalized accuracy
% for i = 1:np
%     prob = eval(problems{i});
%     fprintf(fp, '%s,', prob.name);
%     
%     j = 1;
%     [~, id] = sort(method_results1{j}(i, :));
%     nx = id(mseed);
%     fprintf(fp, '%0.4f & , ',  method_results1{j}(i, nx));
%     fprintf(fp, '%0.4f & , ',  method_results2{j}(i, nx));
%     
%     j = 3;
%     [~, id] = sort(method_results1{j}(i, :));
%     nx = id(mseed);
%     fprintf(fp, '%0.4f & , ',  method_results1{j}(i, nx));
%     fprintf(fp, '%0.4f & , ',  method_results2{j}(i, nx));
%     
%     % conduct significant test cu
%     % 1 means  the former is smaller,  second is bigger,
%     % for accuracy comparison smaller the better
%     [p1,h1,stats1] = ranksum( method_results1{1}(i, :), method_results1{3}(i, :),  'alpha', 0.05, 'tail', 'left');
%     [p2,h2,stats2] = ranksum(method_results1{3}(i, :), method_results1{1}(i, :),  'alpha', 0.05, 'tail', 'left');
%     if h1 == 1 && h2 == 0
%         fprintf(fp, '$\\uparrow$,');
%     elseif h2 == 1 && h1 == 0
%         fprintf(fp, '$\\downarrow$,');
%     else
%         fprintf(fp, '$\\approx$,');
%     end
%     
%     fprintf(fp, '\n');
% end
% fclose(fp);
% 
% %--------------------------------
% 
% 
% filename = strcat('median_accuracy_sigTest_TransferPlus2random_nlvar_', num2str(prob.n_lvar),'.csv');
% fp = fopen(filename, 'w');
% fprintf(fp, 'problems,  transfer to local search+(Upper), (Lower), Random (Upper), (Lower), Significance test \n');
% % normalized accuracy
% for i = 1:np
%     prob = eval(problems{i});
%     fprintf(fp, '%s,', prob.name);
%     
%     j = 2;
%     [~, id] = sort(method_results1{j}(i, :));
%     nx = id(mseed);
%     fprintf(fp, '%0.4f & , ',  method_results1{j}(i, nx));
%     fprintf(fp, '%0.4f & , ',  method_results2{j}(i, nx));
%     
%     j = 3;
%     [~, id] = sort(method_results1{j}(i, :));
%     nx = id(mseed);
%     fprintf(fp, '%0.4f & , ',  method_results1{j}(i, nx));
%     fprintf(fp, '%0.4f & , ',  method_results2{j}(i, nx));
%     
%     % conduct significant test cu
%     % 1 means  the former is smaller,  second is bigger,
%     % for accuracy comparison smaller the better
%     [p1,h1,stats1] = ranksum( method_results1{2}(i, :), method_results1{3}(i, :),  'alpha', 0.05, 'tail', 'left');
%     [p2,h2,stats2] = ranksum(method_results1{2}(i, :), method_results1{3}(i, :),  'alpha', 0.05, 'tail', 'left');
%     if h1 == 1 && h2 == 0
%         fprintf(fp, '$\\uparrow$,');
%     elseif h2 == 1 && h1 == 0
%         fprintf(fp, '$\\downarrow$,');
%     else
%         fprintf(fp, '$\\approx$,');
%     end
%     
%     fprintf(fp, '\n');
% end
% fclose(fp);
% 
% %--------------------------------


norm_up = [];
norm_down = [];
for i = 1:nm
   norm_up = [norm_up; method_results1{i}];
   norm_down = [norm_down; method_results2{i}];
end
norm_up = (norm_up - min(norm_up(:)))/(max(norm_up(:)) - min(norm_up(:)));
norm_down = (norm_down - min(norm_down(:)))/(max(norm_down(:)) - min(norm_down(:)));
norm_accuracy = norm_up + norm_down;

method_normResults = {};

for i = 1:nm
    method_normResults{i} = norm_accuracy((i-1)*np + 1: i*np, :);
end

filename = strcat('median_accuracy_Normlatex_', num2str(prob.n_lvar),'.csv');

fp = fopen(filename, 'w');
fprintf(fp, 'problems,  baseline, transfer to local search, transfer to local search +, decision making,  compare-0, compare-1\n');

for i = 1:np
    prob = eval(problems{i});
    fprintf(fp, '%s &,', prob.name);
    
    for j = 1:nm
        
        [~, id] = sort(method_normResults{j}(i, :));
        nx = id(mseed);
        fprintf(fp, '%0.4f & , ', method_normResults{j}(i, nx));
       
        if j == sigTestIndex %
            % 1 means  the former is smaller,  second is bigger/better, 
            % for hv comparison larger the better
            [p1,h1,stats1] = ranksum( method_results1{j}(i, :), method_results1{j-2}(i, :),  'alpha', 0.05, 'tail', 'left');      
            [p2,h2,stats2] = ranksum( method_results1{j-2}(i, :), method_results1{j}(i, :),  'alpha', 0.05, 'tail', 'left');    
            if h1 == 1 && h2 == 0               
                fprintf(fp, '$\\uparrow_{3}$,');               
            elseif h2==1 && h1 == 0
                fprintf(fp, '$\\downarrow_{3}$,');
            else
                fprintf(fp, '$\\approx_{3}$,');
            end
            
             % 1 means  the former is smaller,  second is bigger/better, 
            % for hv comparison larger the better
            [p1,h1,stats1] = ranksum( method_results1{j}(i, :), method_results1{j-1}(i, :),  'alpha', 0.05, 'tail', 'left');      
            [p2,h2,stats2] = ranksum( method_results1{j-1}(i, :), method_results1{j}(i, :),  'alpha', 0.05, 'tail', 'left');    
            if h1 == 1 && h2 == 0               
                fprintf(fp, '$\\uparrow_{4}$');               
            elseif h2==1 && h1 == 0
                fprintf(fp, '$\\downarrow_{4}$');
            else
                fprintf(fp, '$\\approx_{4}$');
            end
        end
    end
    fprintf(fp, '\\\\\n');

end
fclose(fp);



end

function[] = accuracy_extractionExtension(problems, method, resultfolder, np, seed, mseed, sigTestIndex)

nm = length(method);
method_results1 = {}; % upper save
method_results2 = {}; % lower save


for m = 1:nm
    accuracy_up  = zeros(np, seed);
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
    
    method_results1{m} = accuracy_up;
    method_results2{m} = accuracy_low;
end

filename = strcat('median_accuracy_latex_nlvar_', num2str(prob.n_lvar),'.csv');
fp = fopen(filename, 'w');
fprintf(fp, 'problems,  baseline, _, transfer to local search, _, transfer to local search +, _, decision making, _,  compare-0, compare-1, correction_function, compare-5 \n');

for i = 1:np
    prob = eval(problems{i});
    fprintf(fp, '%s &,', prob.name);
    
    for j = 1:nm
        
        [~, id] = sort(method_results1{j}(i, :));
        nx = id(mseed);
        fprintf(fp, '%0.4f & , ',  method_results1{j}(i, nx));
        fprintf(fp, '%0.4f & , ',  method_results2{j}(i, nx));
 
        if j == sigTestIndex  % accuracy smaller the better
            % 1 means  the former is smaller,  second is bigger/better, 
            % for hv comparison larger the better
            [p1,h1,stats1] = ranksum( method_results1{j}(i, :), method_results1{j-2}(i, :),  'alpha', 0.05, 'tail', 'left');      
            [p2,h2,stats2] = ranksum( method_results1{j-2}(i, :), method_results1{j}(i, :),  'alpha', 0.05, 'tail', 'left');    
            if h1 == 1 && h2 == 0               
                fprintf(fp, '$\\uparrow_{3}$ & , ');               
            elseif h2==1 && h1 == 0
                fprintf(fp, '$\\downarrow_{3}$ &, ');
            else
                fprintf(fp, '$\\approx_{3}$ &, ');
            end
            
            % 1 means  the former is smaller,  second is bigger/better, 
            % for hv comparison larger the better
            [p1,h1,stats1] = ranksum( method_results1{j}(i, :), method_results1{j-1}(i, :),  'alpha', 0.05, 'tail', 'left');      
            [p2,h2,stats2] = ranksum( method_results1{j-1}(i, :), method_results1{j}(i, :),  'alpha', 0.05, 'tail', 'left');    
            if h1 == 1 && h2 == 0               
                fprintf(fp, '$\\uparrow_{4}$& ');               
            elseif h2==1 && h1 == 0
                fprintf(fp, '$\\downarrow_{4}$&');
            else
                fprintf(fp, '$\\approx_{4}$ &');
            end
        end
        
        if j == sigTestIndex+1
            % 1 means  the former is smaller,  second is bigger/better, 
            % accuracy, smaller is better
            [p1,h1,stats1] = ranksum( method_results1{j}(i, :), method_results1{j-1}(i, :),  'alpha', 0.05, 'tail', 'left');      
            [p2,h2,stats2] = ranksum( method_results1{j-1}(i, :), method_results1{j}(i, :),  'alpha', 0.05, 'tail', 'left');    
            if h1 == 1 && h2 == 0               
                fprintf(fp, '$\\uparrow_{5}$ , ');               
            elseif h2 == 1 && h1 == 0
                fprintf(fp, '$\\downarrow_{5}$, ');
            else
                fprintf(fp, '$\\approx_{5}$, ');
            end
        end
    end
    fprintf(fp, '\\\\\n');

end
fclose(fp);


% normalized accuracy 

norm_up = [];
norm_down = [];
for i = 1:nm
   norm_up = [norm_up; method_results1{i}];
   norm_down = [norm_down; method_results2{i}];
end
norm_up = (norm_up - min(norm_up(:)))/(max(norm_up(:)) - min(norm_up(:)));
norm_down = (norm_down - min(norm_down(:)))/(max(norm_down(:)) - min(norm_down(:)));
norm_accuracy = norm_up + norm_down;

method_normResults = {};

for i = 1:nm
    method_normResults{i} = norm_accuracy((i-1)*np + 1: i*np, :);
end

filename = strcat('median_accuracy_Normlatex_', num2str(prob.n_lvar),'.csv');

fp = fopen(filename, 'w');
fprintf(fp, 'problems,  baseline, transfer to local search, transfer to local search +, decision making,  compare-0, compare-1, correctionFunction, compare to 5\n');

for i = 1:np
    prob = eval(problems{i});
    fprintf(fp, '%s &,', prob.name);
    
    for j = 1:nm
        
        [~, id] = sort(method_normResults{j}(i, :));
        nx = id(mseed);
        fprintf(fp, '%0.4f & , ', method_normResults{j}(i, nx));
       
        if j == sigTestIndex %
            % 1 means  the former is smaller,  second is bigger/better, 
            % for hv comparison larger the better
            [p1,h1,stats1] = ranksum( method_results1{j}(i, :), method_results1{j-2}(i, :),  'alpha', 0.05, 'tail', 'left');      
            [p2,h2,stats2] = ranksum( method_results1{j-2}(i, :), method_results1{j}(i, :),  'alpha', 0.05, 'tail', 'left');    
            if h1 == 1 && h2 == 0               
                fprintf(fp, '$\\uparrow_{3}$&,');               
            elseif h2==1 && h1 == 0
                fprintf(fp, '$\\downarrow_{3}$&,');
            else
                fprintf(fp, '$\\approx_{3}$&,');
            end
            
             % 1 means  the former is smaller,  second is bigger/better, 
            % for hv comparison larger the better
            [p1,h1,stats1] = ranksum( method_results1{j}(i, :), method_results1{j-1}(i, :),  'alpha', 0.05, 'tail', 'left');      
            [p2,h2,stats2] = ranksum( method_results1{j-1}(i, :), method_results1{j}(i, :),  'alpha', 0.05, 'tail', 'left');    
            if h1 == 1 && h2 == 0               
                fprintf(fp, '$\\uparrow_{4}$&,');               
            elseif h2==1 && h1 == 0
                fprintf(fp, '$\\downarrow_{4}$&,');
            else
                fprintf(fp, '$\\approx_{4}$&,');
            end
        end
        
        if j == sigTestIndex+1
            % 1 means  the former is smaller,  second is bigger/better, 
            % accuracy, smaller is better
            [p1,h1,stats1] = ranksum( method_results1{j}(i, :), method_results1{j-1}(i, :),  'alpha', 0.05, 'tail', 'left');      
            [p2,h2,stats2] = ranksum( method_results1{j-1}(i, :), method_results1{j}(i, :),  'alpha', 0.05, 'tail', 'left');    
            if h1 == 1 && h2 == 0               
                fprintf(fp, '$\\uparrow_{5}$ , ');               
            elseif h2 == 1 && h1 == 0
                fprintf(fp, '$\\downarrow_{5}$, ');
            else
                fprintf(fp, '$\\approx_{5}$, ');
            end
        end
    end
    fprintf(fp, '\\\\\n');

end
fclose(fp);


end
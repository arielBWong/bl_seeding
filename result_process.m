% generate accuracy analysis

problems = { 'smd1mp(1, 2, 1)' , 'smd2mp(1, 2, 1)',  'smd3mp(1, 2, 1)',  'smd4mp(1, 2, 1)', ....
    'smd5mp(1, 2, 1)' , 'smd6mp(1, 0, 2, 1)', 'smd7mp(1, 2, 1)',  'smd8mp(1, 2, 1)'};


% problems = { 'smd1mp(1, 1, 1)' , 'smd2mp(1, 1, 1)',  'smd3mp(1, 1, 1)',  'smd4mp(1, 1, 1)', ....
%     'smd5mp(1, 1, 1)' , 'smd6mp(1, 0, 1, 1)', 'smd7mp(1, 1, 1)',  'smd8mp(1, 1, 1)'};


method = {'', '0', '1'};
resultfolder = fullfile(pwd, 'resultfolder_2');

np = length(problems);
seed = 29;
mseed = 15;
% [problem, accuracy]


nm = length(method);
method_results1 = {};
method_results2 = {};

FE_results2 = {};

for m = 1:nm
    accuracy_up  = zeros(np, seed);
    accuracy_low = zeros(np, seed);
    lowerFE = zeros(np, seed);
    
    for p = 1:np
        prob = eval(problems{p});
        for s = 1: seed
            filename = strcat('final_accuracy_seed_', num2str(s), '.csv');
            foldername =  strcat(prob.name, method{m});
            savename = fullfile(resultfolder, foldername, filename);
            accuracy = csvread(savename);
            
            accuracy_up(p, s) = accuracy(1);
            accuracy_low(p, s) = accuracy(2);
            
%             filename = strcat('loweFE_seed_', num2str(s), '.csv');
%             foldername =  strcat(prob.name, method{m});
%             savename = fullfile(resultfolder, foldername, filename)
%            lowerFE(p, s) = csvread(savename);
        end
    end
    
    method_results1{m} = accuracy_up;
    method_results2{m} = accuracy_low;
    FE_results2{m} = lowerFE;
end




filename = 'median_accuracy_2.csv';
fp = fopen(filename, 'w');
fprintf(fp, 'problems,  global-, _, median, compare-0, compare-1, seed0, _, median, seed1, _, median, \n');

for i = 1:np
    prob = eval(problems{i});
    fprintf(fp, '%s,', prob.name);
    
    for j = 1:nm
        
        [~, id] = sort(method_results1{j}(i, :));
        nx = id(mseed);
        fprintf(fp, '%0.4f, ', method_results1{j}(i, nx));
        fprintf(fp, '%0.4f,',  method_results2{j}(i, nx));
        fprintf(fp, '%d,',  nx);
        
        
        if j == 1 %
            % 1 means  the former is smaller,  second is bigger/better, 
            % for hv comparison larger the better
            [p1,h1,stats1] = ranksum( method_results1{j}(i, :), method_results1{j+1}(i, :),  'alpha', 0.05, 'tail', 'left');      
            [p2,h2,stats2] = ranksum( method_results1{j+1}(i, :), method_results1{j}(i, :),  'alpha', 0.05, 'tail', 'left');    
            if h1 == 1 && h2 == 0               
                fprintf(fp, 'better-0, ');               
            elseif h2==1 && h1 == 0
                fprintf(fp, 'worse-0, ');
            else
                fprintf(fp, 'same-0, ');
            end
            
             % 1 means  the former is smaller,  second is bigger/better, 
            % for hv comparison larger the better
            [p1,h1,stats1] = ranksum( method_results1{j}(i, :), method_results1{j+2}(i, :),  'alpha', 0.05, 'tail', 'left');      
            [p2,h2,stats2] = ranksum( method_results1{j+2}(i, :), method_results1{j}(i, :),  'alpha', 0.05, 'tail', 'left');    
            if h1 == 1 && h2 == 0               
                fprintf(fp, 'better-1, ');               
            elseif h2==1 && h1 == 0
                fprintf(fp, 'worse-1, ');
            else
                fprintf(fp, 'same-1, ');
            end
        end
    end
    fprintf(fp, '\n');

end
fclose(fp);

% % 
% filename = 'median_FE_3.csv';
% fp = fopen(filename, 'w');
% fprintf(fp, 'problems,  global-, seed0, seed1\n');
% 
% for i = 1:np
%     prob = eval(problems{i});
%     fprintf(fp, '%s,', prob.name);
%     
%     for j = 1:nm
%         
%         [~, id] = sort(FE_results2{j}(i, :));
%         nx = id(mseed);
%         fprintf(fp, '%0.4f, ', FE_results2{j}(i, nx));
%        
%     end
%     fprintf(fp, '\n');
% 
% end
% fclose(fp);


% study the distance metric
distance_metric = {};
distanceN_metric = {};
for i = 1:np
    prob = eval(problems{i});
    filename = strcat(prob.name,'distance_check.csv');
    filename = fullfile('output', filename);
    fp = fopen(filename, 'w');
    fprintf(fp, 'seed,  all_mean, all_median, all_std, all_max, 5_mean, 5_median, 5_std, 5_max\n');
    c = 0;
    for s = 1:seed
        fprintf(fp, '%d,', s);
        filename   = strcat('distanceCheck_seed_', num2str(s), '.csv');
        foldername = strcat(prob.name, method{1});
        savename   = fullfile(resultfolder, foldername, filename);
        d   = csvread(savename);
        da = max(d);
        fprintf(fp, '%0.4f,', mean(d));
        fprintf(fp, '%0.4f,', median(d));
        fprintf(fp, '%0.4f,', std(d));
        fprintf(fp, '%0.4f,', max(d));
        
        filename   = strcat('distanceCheckN_seed_', num2str(s), '.csv');
        foldername = strcat(prob.name, method{1});
        savename   = fullfile(resultfolder, foldername, filename);
        d   = csvread(savename);
        d5 = max(d);
        fprintf(fp, '%0.4f,', mean(d));
        fprintf(fp, '%0.4f,', median(d));
        fprintf(fp, '%0.4f,', std(d));
        fprintf(fp, '%0.4f\n', max(d));
        
        if  d5 <= da
            c = c + 1;
        end
  
    end
    fprintf(fp, '%0.4f', c/seed);
    fclose(fp);
    
end


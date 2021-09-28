% generate accuracy analysis

problems = {'smd5mp(1, 1, 1)' , 'smd7mp(1, 1, 1)',  'smd8mp(1, 1, 1)'};
method = {'', '0', '1'};
resultfolder = fullfile(pwd, 'resultfolder');

np = length(problems);
seed = 11;
m = 6;
% [problem, accuracy]


nm = length(method);
method_results1 = {};
method_results2 = {};

for m = 1:nm
    accuracy_up  = zeros(np, seed);
    accuracy_low = zeros(np, seed);
    
    for p = 1:np
        prob = eval(problems{p});
        for s = 1: seed
            filename = strcat('final_accuracy_seed_', num2str(seed), '.csv');
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




filename = 'median_accuracy.csv';
fp = fopen(filename, 'w');
fprintf(fp, 'problems, infill, _, seed0, _,seed1, _ \n');

for i = 1:np
    prob = eval(problems{i});
    fprintf(fp, '%s,', prob.name);
    
    for j = 1:nm
        
        [~, id] = sort(method_results1{j}(i, :));
        nm = id(m);
        fprintf(fp, '%0.4f, ', method_results1{j}(i, nm));
        fprintf(fp, '%0.4f,',  method_results2{j}(i, nm));
    end
    fprintf(fp, '\n');

end
fclose(fp);



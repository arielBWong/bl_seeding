% move results around 
problem_sets = cell(1, 2);

problem_sets{1} =  {'smd5mp(1, 1, 1)' , 'smd7mp(1, 1, 1)',  'smd8mp(1, 1, 1)', ...
    'smd1mp(1, 1, 1)' , 'smd2mp(1, 1, 1)',  'smd3mp(1, 1, 1)', ...
    'smd4mp(1, 1, 1)', 'smd6mp(1, 0, 1, 1)', };

problem_sets{2} =  {'smd5mp(1, 2, 1)' , 'smd7mp(1, 2, 1)',  'smd8mp(1, 2, 1)', ...
     'smd1mp(1,  2 , 1)' , 'smd2mp(1, 2, 1)',  'smd3mp(1, 2, 1)', ...
     'smd4mp(1, 2, 1)', 'smd6mp(1, 0, 2, 1)', };

 other_threshold = [80,  90];
 
 for is = 1:2
    methods = {'_baseline_ea',  '_seeding_strategy_1'};
    problems = problem_sets{is};
    prob_test = eval(problems{1});
    nv = prob_test.n_lvar;

    np = length(problems);
    seed = 21;
   
    foldername = strcat('resultfolder_trueEval', num2str(nv), '_thr_70');   
    source_resultfolder = fullfile(pwd, foldername);
    
    for it = 1: length(other_threshold)
        foldername = strcat('resultfolder_trueEval', num2str(nv), '_thr_', num2str(other_threshold(it)));
        target_resultfolder = fullfile(pwd, foldername);
        
        % copy source folder and everything under it to target folder
        for ip = 1: length(problems)
            prob = eval(problems{ip});
            for im = 1:length(methods)
                subfolder_name = strcat(prob.name, methods{im});
                target_subfolder_name = fullfile(target_resultfolder,  subfolder_name);
                
                source_subfolder_name = fullfile(source_resultfolder, subfolder_name);               
               
                if ~exist(target_subfolder_name, 'dir')
                    mkdir(target_subfolder_name);
                end
               
                copyfile(source_subfolder_name, target_subfolder_name);
                
                
            end
            
        end
        
    end
    
    
 end
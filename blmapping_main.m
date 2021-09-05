function blmapping_main(problem_str, seed, extended, varargin)
%  This function apply kriging to bl search and bl mapping
%

% parse input
p = inputParser;
addRequired(p,   'problem_str');
addRequired(p,   'seed');
addRequired(p, 'extended');
addParameter(p, 'local_str', 'vanilla');
addParameter(p, 'restart_num', 0);
addParameter(p, 'use_seeding', false);
parse(p, problem_str, seed, extended,  varargin{:});

%--- assign input parameters
prob = p.Results.problem_str;
seed = p.Results.seed;
extended = p.Results.extended;
method_str = p.Results.local_str; % useless 
use_seeding = p.Results.use_seeding;
rx = p.Results.restart_num; 
%---

visualize = false;
if visualize
    fighn = figure('Position', [100 100 800 800]);
    h1     = subplot(2, 2, 1);
    h2     = subplot(2, 2, 3);
    h3     = subplot(2, 2, 4);
end

rng(seed, 'twister');
eval('parameter_load');                                       % IMPORTANT,  variables are load externally
infill_metodstr = 'Believer_nextUpdate';
n_FE = 0;                                                             % record lowerlevel FEs
prob = eval(prob);
uppernext = str2func(infill_metodstr);                % Infill algorithm handle for lower level
method = strcat('local', method_str);
localsearch = false;


%plot_matching(prob);
% lower level settings
llmatch_p.prob                 = prob;
llmatch_p.egostr              = infill_metodstr;
llmatch_p.egofnormstr    = 'normalization_z';    %  Only for place holder/ original Believer and EGO will need
llmatch_p.seed                 = seed;
llmatch_p.localsearch      = localsearch;
llmatch_p.method           = method;
llmatch_p.localmethod    = [];

% pf = prob.upper_pf(100);
% plot_MOupper(prob);
%------------------Process starts--------------------
%---xu initialization
xu = lhsdesign(inisize_u, prob.n_uvar, 'criterion', 'maximin', 'iterations', 1000);
xu = repmat(prob.xu_bl, inisize_u, 1) + repmat((prob.xu_bu - prob.xu_bl), inisize_u, 1) .* xu;

xl	             = [];
low_neval  = 0;

% for testing lower level correlation
train_xl        = lhsdesign(inisize_l, prob.n_lvar, 'criterion','maximin','iterations',1000);
train_xl        = repmat(prob.xl_bl, inisize_l, 1) ...
                        + repmat((prob.xl_bu - prob.xl_bl), inisize_l, 1) .* train_xl;
lower_archive.init_xl = train_xl;
lower_archive.init_fl = [];

%---xu match its xl and evaluate fu
for i = 1:inisize_u
    fprintf('Initialition xu matching process iteration %d\n', i);
    [xl_single, n, ~, lower_archive]   = llmatch_keepdistance(xu(i, :), llmatch_p, 'visualization', false, 'lower_archive', lower_archive);
        
%     lower_archive      = dummy_function(prob, xu(i, :), lower_archive);     
%     xl_single              = prob.get_xlprime(xu(i, :));
%      n                         = 0;
     xl                         = [xl; xl_single];
     n_FE                    = n_FE + n;
end


%---xu evaluation
[fu, cu]                 = prob.evaluate_u(xu, xl);
[fl,  cl]                  = prob.evaluate_l(xu, xl);

% -- archive all evaluated
archive.xu     = xu;
archive.xl       = xl;
archive.fu      = fu;
archive.fl      = fl;
archive.cu     = cu;
archive.cl      = cl;
archive.flag  = ones(inisize_u, 1);

krg_param.GPR_type              = 2;
krg_param.no_trials                = 1;

%---- prepare mapping---
param_ea.popsize = num_pop;
param_ea.gen       = num_gen;
% 

% obj = VideoWriter('moving.avi');
% obj.Quality= 100;
% obj.FrameRate = 25;
% open(obj);

% plot P1
for iter = 1:numiter_u
    fprintf('iterate %d', iter);
    if extended
        [model_xl, arc_trgxl] = creating_blmodel(archive, krg_param);
    else
        model_xl = [];
        arc_trgxl = [];
    end
    [newxu, ~] = Believer_nextExtended(archive, prob, param_ea, krg_param, model_xl, arc_trgxl, extended);
    disp(newxu);
    
    dist  = pdist2(newxu, archive.xu);
    [~, idx] = sort(dist);
    seed_xl = archive.xl(idx(1), :);
    
    if use_seeding
        [newxl, n, ~, lower_archive]   = llmatch_keepdistance(newxu, llmatch_p, 'visualization', false,'seed_xl', seed_xl, ...
            'lower_archive', lower_archive, 'seeding_only', true, 'restartn', rx);   

         if visualize
            seed_xu = archive.xu(idx(1), :);
            seedinglower_plot2D(h1, h2, h3, prob, newxu, newxl, seed_xu, seed_xl, archive, lower_archive, idx)
            
%             for t = 1:30
%                 ff = getframe(gcf);
%                 writeVideo(obj, ff);
%             end    
        end
        
    else
        [newxl, n, ~, lower_archive]   = llmatch_keepdistance(newxu, llmatch_p, 'visualization', false, 'seed_xl', seed_xl, ...
            'lower_archive', lower_archive, 'archive', archive);
        
        % lower_archive   = dummy_function(prob, newxu, lower_archive);
        % newxl = prob.get_xlprime(newxu);
        
    end
    
    low_neval = low_neval + n;
    % newxl = prob.get_xlprime(newxu);
    
    [newfu, newcu] = prob.evaluate_u(newxu, newxl);
    [newfl, newcl] =   prob.evaluate_l(newxu, newxl);
    
  
    archive.xu =  [archive.xu; newxu];
    archive.xl =  [archive.xl; newxl];
    archive.fu =  [archive.fu; newfu];
    archive.fl  =  [archive.fl; newfl];
    archive.cu =  [archive.cu; newcu];
    archive.cl =  [archive.cl; newcl];
    
end
%obj.close();
save_forextended(archive, prob, low_neval, seed, use_seeding, rx);
end

function lower_archive = dummy_function(prob, xu, lower_archive)
n                   = size(lower_archive.init_xl, 1);
xu_init           = repmat(xu, n, 1);
[fl, fc]            = prob.evaluate_l(xu_init, lower_archive.init_xl);
fl                    = fl';
lower_archive.init_fl =[lower_archive.init_fl; fl];
end

%
% for iter = 1:  numiter_u
%     fprintf('upper iter %d\n', iter);
%     if sum(arc_flag) >  size(arc_xu, 1)/2
%
%         if ul_flag < sum(arc_flag)  % rebuild mapping
%             % ul_flag is used to constrol when to update mapping
%             [mapping_bl, mapping_upper, mapping_lower] = create_mappingmodels(arc_xu, arc_xl, arc_fu, arc_fl, arc_flag);
%         end
%
%         % apply mapping based infill
%         [arc_xu, arc_fu, arc_xl, arc_fl, arc_cu, arc_cl, arc_flag] = mapbased_next(archive,...
%             mapping_bl, mapping_upper, mapping_lower,  mapping_lowvar, prob, num_pop, num_gen, norm_str);
%
%
%         % update the monitor to recreate mapping
%         ul_flag = sum(arc_flag);
%         % no evaluation on the lower level              n_FE                     = n_FE + 0;
%
%     else
%         [newxu, ~]            = uppernext(xu, fu,  prob.xu_bu, prob.xu_bl, num_pop, num_gen, cu, norm_str);
%
%         [newxl, n, ~]         = llmatch_keepdistance(newxu, llmatch_p, visual, true);
%         n_FE                     = n_FE + n;
%
%         % (5) flag 0
%         arc_flag = [arc_flag; 1];  % update objective evaluation from lower search
%         % --- true evaluation
%         [newfu, newcu]                        = prob.evaluate_u(newxu, newxl);
%         [newfl, newcl]                          = prob.evaluate_l(newxu, newxl); % redundant but necessary, due to lower level did not return f
%         [arc_xu, arc_xl, arc_fu, arc_fl, arc_cu, arc_cl]                  = add_entry_extened(arc_xu, arc_xl, arc_fu, arc_fl, arc_cu, arc_cl,  newxu, newxl, newfu, newfl, newcu, newcl);
%         [arc_xu, arc_xl, arc_fu, arc_fl, arc_cu, arc_cl, arc_flag]    = keepdistance_upperExtended(arc_xu, arc_xl, arc_fu, arc_fl, arc_cu, arc_cl, arc_flag, prob.xu_bu, prob.xu_bl);
%
%     end
%
%     % visualise
%     if visualize
%         [~, ~, ~, ~, index] = localsolver_startselection(arc_xu,  arc_fu, arc_cu);
%         accuracy_fu = abs( arc_fu(index) -  prob.uopt);
%         plot(iter, accuracy_fu, 'bo');
%         xlim([1,  numiter_u]);
%         hold on;
%
%     end
%
% end

% save_postprocess(arc_xu, arc_xl, arc_fu, arc_fl, arc_cu, arc_cl, arc_flag, prob, seed);


function  save_forextended(archive, prob, low_neval, seed, use_seeding, rx)

xu = archive.xu;
xl = archive.xl;
fu = archive.fu;
fl = archive.fl;
cu = archive.cu;
cl = archive.cl;

resultfolder = fullfile(pwd, 'resultfolder' );
n = exist(resultfolder);
if n ~= 7
    mkdir(resultfolder)
end

if use_seeding
    foldername = strcat(prob.name, '_seeding');
    foldername = strcat(prob.name, num2str(rx));
else
    foldername = prob.name;
end
resultfolder = fullfile(resultfolder, foldername);
n = exist(resultfolder);
if n ~= 7
    mkdir(resultfolder)
end

filename = strcat('xu_seed_', num2str(seed), '.csv');
savename = fullfile(resultfolder, filename);
csvwrite(savename, xu);

filename = strcat('xl_seed_', num2str(seed), '.csv');
savename = fullfile(resultfolder, filename);
csvwrite(savename, xl);

filename = strcat('fu_seed_', num2str(seed), '.csv');
savename = fullfile(resultfolder, filename);
csvwrite(savename, fu);

filename = strcat('fl_seed_', num2str(seed), '.csv');
savename = fullfile(resultfolder, filename);
csvwrite(savename, fl);

filename = strcat('cu_seed_', num2str(seed), '.csv');
savename = fullfile(resultfolder, filename);
csvwrite(savename, cu);

filename = strcat('cl_seed_', num2str(seed), '.csv');
savename = fullfile(resultfolder, filename);
csvwrite(savename, cl);

filename = strcat('loweFE_seed_', num2str(seed), '.csv');
savename = fullfile(resultfolder, filename);
csvwrite(savename, low_neval);

if size(fu, 2) > 1
    % fu nd front
    ndf_indx = Paretoset(fu);
    ndf = fu(ndf_indx, :);
    
    filename = strcat('final_ndfu_seed_', num2str(seed), '.csv');
    savename = fullfile(resultfolder, filename);
    csvwrite(savename, ndf);
    
    
    % fu hv
    pf = prob.upper_pf(100);
    max_f = max(pf, [], 1);
    min_f = min(pf, [], 1);
    norm_fu = (fu - min_f)./(max_f - min_f);
    ref = ones(1, size(pf, 2)) * 1.1;
    
    hv = Hypervolume(norm_fu, ref);
    filename = strcat('final_ndhv_seed_', num2str(seed), '.csv');
    savename = fullfile(resultfolder, filename);
    csvwrite(savename, hv);
    
    % fu hv
    igd = mean(min(pdist2(pf, ndf), [], 2));
    filename = strcat('final_ndigd_seed_', num2str(seed), '.csv');
    savename = fullfile(resultfolder, filename);
    csvwrite(savename, igd);
else
    [best_x, best_f, best_c, s, index] =  localsolver_startselection(xu,  fu, cu);
    ulp = abs(prob.fu_prime - best_f);
    llp = abs(prob.fl_prime - fl(index) );
    filename = strcat('final_accuracy_seed_', num2str(seed), '.csv');
    savename = fullfile(resultfolder, filename);
    csvwrite(savename, [ulp, llp]);
    
    
end


end

function plotlower(prob, newxu)
n = 50;
xl1 = linspace( prob.xl_bl(1), prob.xl_bu(1), n);
xl2 = linspace(prob.xl_bl(2), prob.xl_bu(2), n);

[xl11, xl22] = meshgrid(xl1, xl2);
fl_pred = zeros(n, n);
fl_real = zeros(n, n);
for i = 1: n
    for j =1: n
        xll = [xl11(i, j), xl22(i, j)];
        
        fl_real(i, j) = prob.evaluate_l(newxu, xll);
    end
end

subplot(1, 1, 1);
surf(xl1, xl2, fl_real , 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
pause(1);
close;
end

function plot_matchingaccuracy(prob, arc_xu, arc_xl, arc_fl, arc_flag, mapping_bl, mapping_upper, mapping_lower,  mapping_lowvar)

n = 50;
xu1 = linspace( prob.xu_bl(1), prob.xu_bu(1), n);
xu2 = linspace(prob.xu_bl(2), prob.xu_bu(2), n);

fl_training = arc_fl(arc_flag==1,:);

[xu11, xu22] = meshgrid(xu1, xu2);
fl= zeros(n, n);
fl_mapping= zeros(n, n);
for i = 1: n
    for j =1: n
        xuu = [xu11(i, j), xu22(i, j)];
        
        [xl(i, j), ~] = Predict_DACE(mapping_lowvar{2}, xuu, 0);
        xl_real = prob. get_xlprime(xuu);
    end
end
subplot(2, 2 , 1);
surf(xu11, xu22, xl, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;

subplot(2, 2 , 2);
surf(xu11, xu22, xl_real, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;




%
% subplot(1, 3 , 2);
% surf(xu11, xu22, fl_mapping, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;
%
%
% %
% %
% [xu11, xu22] = meshgrid(xu1, xu2);
% fl_mapping= zeros(n, n);
% for i = 1: n
%     for j =1: n
%         xuu = [xu11(i, j), xu22(i, j)];
%
%         % init for local search
%         % dist = pdist2(xuu, arc_xu);
%         % [~, ind] =sort(dist);
%         % x0 = arc_xl(ind(1), :);
%         x0= [];
%         for k = 1:length(mapping_lowvar)
%             [x, ~] = Predict_DACE(mapping_lowvar{k}, xuu, 0);
%             x0 = [x0, x];
%         end
%         x0 = boundary_check(x0, prob.xl_bu, prob.xl_bl);
%         %
%         %[fl_bound, ~] = Predict_DACE(mapping_bl{1}, xuu, 0);
%
% %         options = optimset('Algorithm','sqp','Display','off');
% %         [xl, ~, ~, out] = fmincon(@(xl)obj(xuu, xl, mapping_upper), ...
% %             x0, ... %x0
% %             [], ... % A
% %             [], ... % b
% %             [], ... % Aeq
% %             [], ... % beq
% %             prob. xl_bl, ...
% %             prob. xl_bu, ...
% %             @(xl)cons(xuu, xl, mapping_lower, fl_bound),...
% %             options);
%
%         fl_mapping(i, j) = prob.evaluate_l(xuu, x0);
%     end
% end
% subplot(1, 3, 3);
% surf(xu11, xu22, fl_mapping, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;
%
%

end


function plot_matching(prob)
n = 50;
xu= linspace( prob.xu_bl(1), prob.xu_bu(1), n);
xl =linspace( prob.xl_bl(1), prob.xl_bu(1), n);

[xu11, xu22] = meshgrid(xu, xl);
fu= zeros(n, n);
for i = 1: n
    for j =1: n
        fu(i, j) =prob.evaluate_u(xu11(i, j), xu22(i, j));
        
    end
end

surf(xu11, xu22, fu, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;

a = 0;

end

function plot_MOupper(prob)
n = 50;
xu1= linspace( prob.xu_bl(1), prob.xu_bu(1), n);
xu2=linspace( prob.xu_bl(2), prob.xu_bu(2), n);

[xu11, xu22] = meshgrid(xu1, xu2);
fu= zeros(n, n);
for i = 1: n
    for j =1: n
        xu = [xu11(i, j), xu22(i, j)];
        xl = prob.get_xlprime(xu);
        f =prob.evaluate_u(xu, xl);
        fu(i, j) = f(1);
    end
end

surf(xu11, xu22, fu, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;
xlabel('xu1');
ylabel('xu2');
a = 0;

end






function f = obj(xu, xl, model_upper)
x = [xu, xl];
[f, ~] = Predict_DACE(model_upper{1}, x, 0);
end

function [c, ceq] = cons(xu, xl, model_lower, fl_bound)
% create a mapping constraints
ceq = [];
x = [xu, xl];
[fl, ~] = Predict_DACE(model_lower{1}, x, 0);
c = fl - fl_bound;
end


function save_postprocess(xu, xl, fu, fl, cu, cl, arc_flag, prob, seed, extended)
resultfolder = fullfile(pwd, 'resultfolder' );
n = exist(resultfolder);
if n ~= 7
    mkdir(resultfolder)
end
if extended == true
    folder = strcat(prob.name, '_ext');
else
    folder = prob.name;
end
resultfolder = fullfile(resultfolder, folder);
n = exist(resultfolder);
if n ~= 7
    mkdir(resultfolder)
end



filename = strcat('xu_seed_', num2str(seed), '.csv');
savename = fullfile(resultfolder, filename);
csvwrite(savename, xu);

filename = strcat('xl_seed_', num2str(seed), '.csv');
savename = fullfile(resultfolder, filename);
csvwrite(savename, xl);

filename = strcat('fu_seed_', num2str(seed), '.csv');
savename = fullfile(resultfolder, filename);
csvwrite(savename, fu);

filename = strcat('fl_seed_', num2str(seed), '.csv');
savename = fullfile(resultfolder, filename);
csvwrite(savename, fl);

filename = strcat('cu_seed_', num2str(seed), '.csv');
savename = fullfile(resultfolder, filename);
csvwrite(savename, cu);


filename = strcat('cl_seed_', num2str(seed), '.csv');
savename = fullfile(resultfolder, filename);
csvwrite(savename, cl);


filename = strcat('flag_seed_', num2str(seed), '.csv');
savename = fullfile(resultfolder, filename);
csvwrite(savename, arc_flag);

[~, ~, ~, ~, index] = localsolver_startselection(xu,  fu, cu);
accuracy_fu = abs( fu(index) -  prob.uopt);
accuracy_fl = abs(fl(index) - prob.lopt);

accuracy = [accuracy_fu, accuracy_fl];
filename = strcat('accuracy_seed_', num2str(seed), '.csv');
savename = fullfile(resultfolder, filename);
csvwrite(savename, accuracy);

end




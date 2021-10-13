function blmapping_main(problem_str, seed, extended, varargin)
%  This function apply kriging to bl search and bl mapping
%

% parse input
p = inputParser;
addRequired(p, 'problem_str');
addRequired(p, 'seed');
addRequired(p, 'extended');
addParameter(p,'local_str', 'vanilla');
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
    fighn  = figure('Position', [100 100 800 800]);
    h1     = subplot(2, 2, 1);
    h2     = subplot(2, 2, 3);
    h3     = subplot(2, 2, 4);
end

rng(seed, 'twister');
eval('parameter_load2');                                       % IMPORTANT,  variables are load externally
infill_metodstr = 'Believer_nextUpdate';
n_FE = 0;                                                             % record lowerlevel FEs
prob = eval(prob);
uppernext = str2func(infill_metodstr);                % Infill algorithm handle for lower level
method = strcat('local', method_str);
localsearch = false;
normhn = str2func('normalization_z');


%plot_matching(prob);
% lower level settings
llmatch_p.prob                = prob;
llmatch_p.egostr              = infill_metodstr;
llmatch_p.egofnormstr         = 'normalization_z';    %  Only for place holder/ original Believer and EGO will need
llmatch_p.seed                = seed;
llmatch_p.localsearch      = localsearch;
llmatch_p.method           = method;
llmatch_p.localmethod    = [];

% pf = prob.upper_pf(100);
% plot_MOupper(prob);
%------------------Process starts--------------------
%---xu initialization
xu = lhsdesign(inisize_u, prob.n_uvar, 'criterion', 'maximin', 'iterations', 1000);
xu = repmat(prob.xu_bl, inisize_u, 1) + repmat((prob.xu_bu - prob.xu_bl), inisize_u, 1) .* xu;

xl	       = [];
low_neval  = 0;

% for testing lower level correlation
train_xl        = lhsdesign(inisize_l, prob.n_lvar, 'criterion','maximin','iterations',1000);
train_xl        = repmat(prob.xl_bl, inisize_l, 1) ...
                        + repmat((prob.xl_bu - prob.xl_bl), inisize_l, 1) .* train_xl;
lower_archive.init_xl = train_xl;
lower_archive.init_fl = [];

%---xu match its xl and evaluate fu
% artifical_flag = [];
for i = 1:inisize_u
    fprintf('Initialition xu matching process iteration %d\n', i);
    [xl_single, n, flag, lower_archive]   = llmatch_keepdistance(xu(i, :), llmatch_p, 'visualization', false, 'lower_archive', lower_archive);
    % artifical_flag     = [artifical_flag; flag];
     % xl_single         = prob.get_xlprime(xu(i, :));
     % n                 = 0;
     xl                  = [xl; xl_single];
     n_FE                = n_FE + n;
end


%---xu evaluation
[fu, cu]                 = prob.evaluate_u(xu, xl);
[fl,  cl]                = prob.evaluate_l(xu, xl);

% -- archive all evaluated
archive.xu     = xu;
archive.xl     = xl;
archive.fu     = fu;
archive.fl     = fl;
archive.cu     = cu;
archive.cl     = cl;

% archive.artifical_flag = artifical_flag;

krg_param.GPR_type = 2;
krg_param.no_trials = 1;

%---- prepare mapping---
param_ea.popsize = num_pop;
param_ea.gen       = num_gen;

map_param.GPR_type = 2;
map_param.no_trials = 1;


% plot P1
for iter = 1:numiter_u
    fprintf('xu infill iterate %d \n', iter);
     
%     mapping_archive.x   = archive.xu;
%     mapping_archive.muf = archive.fl;
%     mapping_archive.mug = archive.cl;
%     mapping_bl = Train_GPR(archive.xu , archive.fl, map_param);
%     
    [newxu, ~] = Believer_nextUpdate(xu, fu, prob.xu_bu, prob.xu_bl, ...
                        num_pop, num_gen, cu, normhn);
                                        
    % [pred_fl, pred_mu] = Predict_GPR(mapping_bl, newxu, map_param, mapping_archive);
    % disp(newxu);

    dist = pdist2(newxu, archive.xu);
    [~, idx] = sort(dist);
    seed_xl = archive.xl(idx(1), :);
    
    if use_seeding
        [newxl, n, ~, lower_archive]   = llmatch_keepdistance(newxu, llmatch_p, 'visualization', false,'seed_xl', seed_xl, ...
            'lower_archive', lower_archive, 'seeding_only', true, 'restartn', rx);      
    else
        [newxl, n, flag, lower_archive]   = llmatch_keepdistance(newxu, llmatch_p, 'visualization', false, 'seed_xl', seed_xl, ...
            'lower_archive', lower_archive, 'archive', archive); %, 'pred_fl',pred_fl, 'pred_mu', pred_mu);
    end
    
   
    low_neval = low_neval + n;
    % newxl = prob.get_xlprime(newxu);
    
    [newfu, newcu] = prob.evaluate_u(newxu, newxl);
    [newfl, newcl] = prob.evaluate_l(newxu, newxl);
    
    if visualize
        map_blreal(h1, prob);
        map_blpred(h2, mapping_bl, prob, map_param, mapping_archive)
        vis_upper(h3, prob, archive, newxu, newfu);
    end
    
    
    % distance control
    distcheck_xu = [xu; newxu];
    distcheck_fu = [fu; newfu];
    distcheck_cu = [cu; newcu];
    
    dist_id = keepdistance_returnID(distcheck_xu, distcheck_fu, distcheck_cu, prob.xu_bu, prob.xu_bl);
    
    % expand
    xu = [xu; newxu]; xl = [xl; newxl]; fu = [fu; newfu];
    fl = [fl; newfl]; cu = [cu; newcu]; cl = [cl; newcl];
    % artifical_flag = [artifical_flag; flag]; %test step
    % eliminate close
    if dist_id > 0
        
        fprintf('upper to delete %d \n', dist_id);
        fprintf('upper total size %d \n',size(xu, 1));
        
        xu(dist_id, :) = []; xl(dist_id, :) = []; fu(dist_id, :) = [];
        fl(dist_id, :) = [];  % artifical_flag(dist_id, :) = [];
        if ~isempty(cu) 
            cu(dist_id, :) = [];
        end
        if ~isempty(cl)
            cl(dist_id, :) = [];
        end
    end
    
   % save 
    archive.xu =  xu; archive.xl =  xl;  archive.fu =  fu;
    archive.fl = fl;  archive.cu =  cu; archive.cl =  cl;
    % archive.artifical_flag = artifical_flag;

    
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




function  save_forextended(archive, prob, low_neval, seed, use_seeding, rx)

xu = archive.xu;
xl = archive.xl;
fu = archive.fu;
fl = archive.fl;
cu = archive.cu;
cl = archive.cl;

name = strcat('resultfolder_' , num2str(prob.n_lvar));
resultfolder = fullfile(pwd, name);
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


% check which one lower level goes wrong
% whether it starts from where local search is misled
fl_pp = [];
for i = 1:size(fl, 1)
    xl_prime = prob.get_xlprime(xu(i, :));   
    fl_prime = prob.evaluate_l(xu(i, :), xl_prime);
    fl_pp = [fl_pp; fl_prime];    
end

fl = [fl, fl_pp, abs(fl - fl_pp)];
filename = strcat('loweAccuracy_seed_', num2str(seed), '.csv');
savename = fullfile(resultfolder, filename);
csvwrite(savename, fl);


fu_pp = [];
for i = 1:size(fl, 1)
    xl_prime = prob.get_xlprime(xu(i, :));   
    fu_prime = prob.evaluate_u(xu(i, :), xl_prime);
    fu_pp = [fu_pp; fu_prime];    
end

fu = [fu, fu_pp, abs(fu - fu_pp)];
filename = strcat('upperAccuracy_seed_', num2str(seed), '.csv');
savename = fullfile(resultfolder, filename);
csvwrite(savename, fl);


end


function map_blreal(h1, prob)
cla(h1);
n = 50;
x1_tst = linspace( prob.xu_bl(1), prob.xu_bu(1), n);
x2_tst = linspace(prob.xu_bl(2), prob.xu_bu(2), n);

[msx1, msx2]        = meshgrid(x1_tst, x2_tst);
msx11 = msx1(:);
msx22 = msx2(:);
xu = [msx11, msx22];
xl = prob.get_xlprime(xu);
fl = prob.evaluate_l(xu, xl);
fl = reshape(fl, [n, n]);

axes(h1);

surf(msx1, msx2, fl, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;


end
function map_blpred(subplothn, mapping_bl, prob, map_param, mapping_archive)
cla(subplothn);
n = 50;
x1_tst = linspace( prob.xu_bl(1), prob.xu_bu(1), n);
x2_tst = linspace(prob.xu_bl(2), prob.xu_bu(2), n);

[msx1, msx2]        = meshgrid(x1_tst, x2_tst);
msx11 = msx1(:);
msx22 = msx2(:);
xu = [msx11, msx22];
fl =  Predict_GPR(mapping_bl, xu, map_param, mapping_archive);
fl = reshape(fl, [n, n]);
axes(subplothn);
surf(msx1, msx2, fl, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;


end


function vis_upper(subplothn, prob, archive, newxu, newfu)
cla(subplothn);
n = 50;
x1_tst = linspace( prob.xu_bl(1), prob.xu_bu(1), n);
x2_tst = linspace(prob.xu_bl(2), prob.xu_bu(2), n);

[msx1, msx2]        = meshgrid(x1_tst, x2_tst);
msx11 = msx1(:);
msx22 = msx2(:);
xu = [msx11, msx22];
xl = prob.get_xlprime(xu);
fu = prob.evaluate_u(xu, xl);
fu = reshape(fu, [n, n]);

axes(subplothn);

surf(msx1, msx2, fu, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;
scatter3(archive.xu(:, 1), archive.xu(:, 2), archive.fu, 80, 'r', 'filled');
scatter3(newxu(:, 1), newxu(:, 2), newfu, 80, 'g', 'filled');



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






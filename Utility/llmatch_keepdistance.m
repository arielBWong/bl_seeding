function[match_xl, n_fev, flag, lower_archive] = llmatch_keepdistance(upper_xu, lower_parameters, varargin)
% method of searching for a match xl for xu.
% Problem(Prob) definition require certain formation for bilevel problems
% evaluation method should be of  form 'evaluation_l(xu, xl)'
%--------------------------------------------------------------------------
% distribute parameters

% parse input
%-------
p = inputParser;
addRequired(p, 'upper_xu');
addRequired(p, 'lower_parameters');
addParameter(p,'visualization', false);
addParameter(p,'seed_xl', []);
addParameter(p,'lower_archive', []);                   %  probing archive
addParameter(p,'seeding_only', false);                            %  conduct seeding controller
addParameter(p,'restartn', 0);                             %  additional restart times
addParameter(p,'archive', []);                             %   external archive
parse(p,upper_xu, lower_parameters, varargin{:});

llmatch_p      = p.Results.lower_parameters;
xu                 = p.Results.upper_xu;
seed_xl         = p.Results.seed_xl;
visualization = p.Results.visualization;
seeding_only= p.Results.seeding_only;
rn                 = p.Results.restartn;
lower_archive=p.Results.lower_archive;
archive          = p.Results.archive;
%--------

prob                   = llmatch_p.prob;
num_pop            = llmatch_p.num_pop;
num_gen             = llmatch_p.num_gen;
propose_nextx    = llmatch_p.egostr;
init_size               = llmatch_p.egoinitsize;
iter_size               = llmatch_p.egoitersize;
norm_str             = llmatch_p.egofnormstr;
localsearch          = llmatch_p.localsearch;
method               = llmatch_p.method;
seed                    = llmatch_p.seed;


l_nvar                  = prob.n_lvar;
upper_bound      = prob.xl_bu;
lower_bound      = prob.xl_bl;
xu_init                 = repmat(xu, init_size, 1);

%-- lower level initialization
if ~isempty(lower_archive)                % for lower level archive 
    train_xl = lower_archive.init_xl;
else
    train_xl        = lhsdesign(init_size,l_nvar,'criterion','maximin','iterations',1000);
    train_xl        = repmat(lower_bound, init_size, 1) ...
                           + repmat((upper_bound - lower_bound), init_size, 1) .* train_xl;
end

% include seeding xl
if ~isempty(seed_xl)  && seeding_only                                 %  only local search with seeding
    if rn>0
        maxFE = 100; else
        maxFE = 200;
    end
    [seed_fl, seed_fc]                = prob.evaluate_l(xu, seed_xl);
    [match_xl, flag, num_eval] = ll_localsearch(seed_xl, seed_fl, seed_fc, true, xu, prob, maxFE);
    n_fev = num_eval + 1;
    [match_fl, match_cl]          = prob.evaluate_l(xu, match_xl); % additional lazy step, can be extracted from local search results
    % lower_archive = varargin{4};
    
    % random restart other local seasrch points
    if rn>0 % number of restart
        restart_xl  = lhsdesign(rn, prob.n_lvar, 'criterion','maximin','iterations',1000);
        restart_xl  = repmat(prob.xl_bl, rn, 1) ...
            + repmat((prob.xl_bu - prob.xl_bl), rn, 1) .* restart_xl;
        
        additional_searchxl = [];
        additional_searchfl = [];
        additional_searchcl = [];
        for i = 1:rn
            [restart_fl, restart_cl] = prob.evaluate_l(xu, restart_xl(i, :));
            [match_xlrestart, flag, num_eval] = ll_localsearch(restart_xl(i, :), restart_fl , restart_cl, true, xu, prob, maxFE);
            n_fev = n_fev + num_eval + 1;
            [match_flrestart, match_clrestart] = prob.evaluate_l(xu, match_xlrestart);
            additional_searchxl = [additional_searchxl; match_xlrestart];
            additional_searchfl = [additional_searchfl; match_flrestart];
            additional_searchcl = [additional_searchcl; match_clrestart];
        end
        
        % select with previous search results with seeding
        additional_searchxl  = [additional_searchxl; match_xl];
        additional_searchfl   = [additional_searchfl; match_fl];
        additional_searchcl  = [additional_searchcl; match_cl];
        
        [match_xl, ~ , ~, ~] =  localsolver_startselection(additional_searchxl, additional_searchfl, additional_searchcl);
    end
    
    % the following two lines are not needed
    % just for checking correlation
    train_xl = lower_archive.init_xl;
    [train_fl, train_fc] = prob.evaluate_l(xu_init, train_xl);
    lower_archive.init_fl =[lower_archive.init_fl; train_fl'];
    return
    
end

% evaluate/get training fl from xu_init and train_xl
% compatible with non-constriant problem
[train_fl, train_fc] = prob.evaluate_l(xu_init, train_xl);
lower_archive.init_fl =[lower_archive.init_fl; train_fl']; 


% decide whether skip infill steps
if ~isempty(lower_archive) && ~isempty(archive)  % when archive is passed in, means can do closeness check
    dist  = pdist2(xu , archive.xu);      % this xu is upper level new infill xu, not added into archive
    [~, idx] = sort(dist);
    r = corr(train_fl, lower_archive.init_fl(idx(1), :)');

    
    if r>0.95 % skip infill
        maxFE = 180;
       [seed_fl, seed_fc]                = prob.evaluate_l(xu, seed_xl);
       [match_xl, ~, num_eval]     = ll_localsearch(seed_xl, seed_fl, seed_fc, true, xu, prob, maxFE);
       n_fev = num_eval + 1;
       [match_fl, match_cl]          = prob.evaluate_l(xu, match_xl);  % additional lazy step, can be extracted from local search results
       additional_searchxl           = [seed_xl; match_xl]; % variable name from copy paste 
       additional_searchfl           = [seed_fl; match_fl];
       additional_searchcl           = [seed_fc; match_cl];
        
        [match_xl, ~ , flag, ~] =  localsolver_startselection(additional_searchxl, additional_searchfl, additional_searchcl);
        n_fev = n_fev + size(train_xl, 1);
        return
    end
end

if false                       % this is for insert seeding into infill not necessary
    train_xl = [train_xl; seed_xl];
    [seed_fl, seed_fc] = prob.evaluate_l(xu, seed_xl);
    train_fl = [train_fl; seed_fl];
    train_fc = [train_fc; seed_fc];
   iter_size = iter_size -1;
   init_size = init_size + 1;
end


% General infill process
arc_xl               = train_xl;
arc_fl               = train_fl;
arc_cl               = train_fc;

% call EIM to expand train xl one by one
nextx_hn             = str2func(propose_nextx);
normhn               = str2func(norm_str);
if visualization
    fighn            = figure(1);
end

if localsearch
    localmodeling    = str2func(llmatch_p.localmethod);
end


iter = 1;
while size(arc_xl, 1) <= iter_size + init_size
    
    if size(arc_xl, 1) == iter_size + init_size
        break;
    end
    % fprintf('iteration %d\n', iter);
    
    % evaluate dace compatibility
    [train_xl, train_fl, train_fc, ~] ...
        = data_prepare(train_xl, train_fl, train_fc);
    % fprintf('surrogate training data size %d\n', size(train_xl, 1));
    % if eim propose next xl
    % lower level is single objective so no normalization method is needed
    [new_xl, infor]  = nextx_hn(train_xl, train_fl, upper_bound, lower_bound, ...
        num_pop, num_gen, train_fc, normhn, prob);
    
    % local search on surrogate
    % evaluate next xl with xu
    [new_fl, new_fc] = prob.evaluate_l(xu, new_xl);
    
    % visualization
    if l_nvar == 1 && visualization
        processplot1d(fighn, train_xl, train_fl, infor.krg, prob, new_xl); end
    if l_nvar == 2 && visualization %&& localsearch == false
        plotlocalKprocess2D(fighn, prob, train_xl, train_fl, infor.krg, infor.krgc, ...
            new_xl, new_fl, infor.arc_obj, infor.arc_c, xu);
    end
    
    %-----------------------------
    train_xl = [train_xl; new_xl];
    train_fl = [train_fl; new_fl];
    train_fc = [train_fc; new_fc];  % compatible with nonconstraint
    
    [train_xl,train_fl, train_fc]...
             = keepdistance(train_xl,train_fl, train_fc, prob.xl_bu, prob.xl_bl);
    
    arc_xl   = [arc_xl; new_xl];
    arc_fl   = [arc_fl; new_fl];
    arc_cl   = [arc_cl; new_fc];

    
    if size(arc_xl, 1) == iter_size + init_size
        break;
    end
    
    
    if mod(iter, 5) == 0 && localsearch
        % extract local points -- search xl -- evaluate fl--
        % update archives
              
        [localxl, localfl, localcl, krg, krgc, arc_obj, arc_c] ...
            = localmodeling(train_xl, train_fl, train_fc, prob);
        
        
        if ~isempty(krg)
            local_ub                    = max(localxl, [], 1);
            local_lb                     = min(localxl, [], 1);
            
            new_localxl              = localsurrgoate_search(krg, krgc, local_ub, local_lb,  arc_obj, arc_c,...
                num_pop, num_gen, prob.xl_bl, prob.xl_bu);
            [new_localfl, new_localcl]  = prob.evaluate_l(xu, new_localxl);
            
            arc_xl                      = [arc_xl; new_localxl];
            arc_fl                      = [arc_fl; new_localfl];
            arc_cl                      = [arc_cl; new_localcl];
            
            if size(arc_xl, 1) == iter_size + init_size
                break;
            end
            
            %-------------------------------------------------------            
            if true %include_local
            train_xl                    = [train_xl; new_localxl];
            train_fl                    = [train_fl; new_localfl];
            train_fc                    = [train_fc; new_localcl];  % compatible with nonconstraint
            
            [train_xl,train_fl, train_fc]...
                                           = keepdistance(train_xl,train_fl, train_fc, prob.xl_bu, prob.xl_bl);
            end
            % -------------------------------------------------------
            
            
            if l_nvar == 1 && visualization
                plotlocalKprocess(fighn, prob, arc_xl, arc_fl, localxl, localfl, krg, ...
                    new_localxl, new_localfl, arc_obj);
            end
            if l_nvar == 2 && visualization
                plotlocalKprocess2D(fighn, prob, localxl, localfl, krg, krgc, ...
                    new_localxl, new_localfl, arc_obj, arc_c, xu, true);
            end
        else
            prob_name = prob.name;
            fprintf('%s, seed %d, iteration %d local search fails on local krg training failed \n', prob_name, seed, iter);
        end
        
    end
    iter = iter + 1;
end
%fprintf('true iteration is %d\n', iter);

[best_x, best_f, best_c, s] =  localsolver_startselection(arc_xl, arc_fl, arc_cl);
nolocalsearch = false;
if nolocalsearch
    match_xl = best_x;
    n_fev    = size(arc_xl, 1);
    flag     = s;
else
    if size(train_fl, 2)> 1
        error('local search does not apply to MO');
    end
    maxFE = 160;
    [match_xl, ~, num_eval] = ll_localsearch(best_x, best_f, best_c, s, xu, prob, maxFE);
    n_global                   = size(train_xl, 1);
    n_fev                      = n_global +num_eval;       % one in a population is evaluated
   
    % --- avoid sqp overshooting problem 
    [match_fl, match_cl]         = prob.evaluate_l(xu, match_xl); % lazy step, no FE should be counted here
    additional_searchxl           = [best_x; match_xl]; % variable name from copy paste
    additional_searchfl           = [best_f;  match_fl];
    additional_searchcl           = [best_c; match_cl];
    
    [match_xl, ~ , flag, ~] =  localsolver_startselection(additional_searchxl, additional_searchfl, additional_searchcl);
    
end




% save lower level
% llcmp = true;
llcmp = false;
if llcmp
    % only for SO
    if size(train_fl, 2) ==  1
        arc_xl                = [arc_xl; match_xl];
        [local_fl, local_fc]  = prob.evaluate_l(xu, match_xl);
        arc_fl                = [arc_fl; local_fl];
        arc_cl                = [arc_cl; local_fc];
    end
    savelower(prob,arc_xl, arc_fl, arc_cl, method, seed);
end
end



function savelower(prob, x, f, c, method, seed)
num = length(prob.xl_bl);
savepath = strcat(pwd, '\resultfolder_ll\');
n = exist(savepath);
if n ~= 7
    mkdir(savepath)
end


savepath = strcat(pwd, '\resultfolder_ll\', prob.name, '_', num2str(num) ,'_',method);
n = exist(savepath);
if n ~= 7
    mkdir(savepath)
end

% check and process for multiple objective
% save nd front only
if size(f, 2)> 1
    % extract nd front
    num_con = size(c, 2);
    if ~isempty(c) % constraint problems
        index_c = sum(c <= 0, 2) == num_con;
        if sum(index_c) ~=0  % exist feasible,
            feasible_y = f(index_c, :);
            feasible_x = x(index_c, :);
            feasible_c = c(index_c, :);
            
            nd_index = Paretoset(feasible_y);
            f_nd = feasible_y(nd_index, :);
            x_nd = feasible_x(nd_index, :);
            c_nd = feasible_c(nd_index, :);
            
        else
            f_nd =[];
            x_nd = [];
            c_nd = [];
        end
    else % non constraint upper problem/ nd always exists
        nd_index = Paretoset(f);
        f_nd = f(nd_index, :);
        x_nd = x(nd_index, :);
        c_nd = [];
    end
    %---------------
    % save nd front
    savename_xu = strcat(savepath, '\xl_', num2str(seed),'.csv');
    savename_fu = strcat(savepath, '\fl_', num2str(seed),'.csv');
    savename_fc = strcat(savepath, '\cl_', num2str(seed),'.csv');
    
    csvwrite(savename_xu, x_nd);
    csvwrite(savename_fu, f_nd);
    csvwrite(savename_fc, c_nd);
    
    return
end

savename_xu = strcat(savepath, '\xl_', num2str(seed),'.csv');
savename_fu = strcat(savepath, '\fl_', num2str(seed),'.csv');
savename_fc = strcat(savepath, '\cl_', num2str(seed),'.csv');

csvwrite(savename_xu, x);
csvwrite(savename_fu, f);
csvwrite(savename_fc, c);
end

function f = denormzscore(trainy, fnorm)
[~, y_mean, y_std] = zscore(trainy);
f = fnorm * y_std + y_mean;
end

% ---demo test: forrestor
function[] = processplot1d(fighn, trainx, trainy, krg, prob, before)

param.GPR_type   = 2;
param.no_trials  = 5;

train_y_norm     = normalization_z(trainy);
arc_obj.x        = trainx;
arc_obj.muf      = train_y_norm;

clf(fighn);
% (1) create test
testdata        = linspace(prob.xl_bl, prob.xl_bu, 10000);
testdata        = testdata';

% (2) predict
[fpred, sig]    = Predict_GPR(krg, testdata, param, arc_obj);

yyaxis left;
fpred           = denormzscore(trainy, fpred);
plot(testdata, fpred, 'r', 'LineWidth', 1); hold on;

% sig variation
y1              = fpred + sig;
y2              = fpred - sig;
y               = [y1', fliplr(y2')];
x               = [testdata', fliplr(testdata')];
fill(x, y, 'r', 'FaceAlpha', 0.1, 'EdgeColor','none'); hold on;

% (3) real
[freal, sig]    = prob.evaluate_l([], testdata);
plot(testdata, freal, 'b', 'LineWidth', 2);hold on;

% (4) scatter train
scatter(trainx, trainy, 80, 'ko', 'LineWidth', 2);

%---
newy           = prob.evaluate_l([], before);
scatter(before, newy, 80, 'ko', 'filled');

% (5) calculate EI and plot
yyaxis right;
[train_ynorm, ~, ~] = zscore(trainy);

ynorm_min           = min(train_ynorm);



fit          = EIM_evaldaceUpdate(testdata, ynorm_min,  arc_obj, krg);
fit          = -fit;
plot(testdata, fit, 'g--');

data         = [trainx, trainy];
save('data', 'data');
pause(0.5);
end

function plotEIMprocess(fighn, trgx, trgf, nextx, nextf, prob, krg)
clf(fighn);
% (1) create test
testdata = linspace(prob.xl_bl, prob.xl_bu, 100);
testdata = testdata';

[freal, ~]= prob.evaluate_l([], testdata);
plot(testdata, freal, 'b', 'LineWidth', 1);hold on;

% (2) scatter training data
scatter(trgx, trgf, 80, 'ko', 'filled');

% (3) scatter test data
scatter(nextx, nextf, 80, 'ro',  'filled');

% (4) scatter krg prediction
[fpred, ~] = surrogate_predict(testdata, krg, true);
fpred = denormzscore(trgf, fpred);
plot(testdata, fpred, 'r--', 'LineWidth', 1);hold on;
pause(1);
end

function plotlocalKprocess(fighn, prob, arc_xl, arc_fl, local_xl, local_fl, krg, ...
    localnewx, localnewf, arc_obj)

clf(fighn);
% (1) create test
testdata    = linspace(prob.xl_bl, prob.xl_bu, 1000);
testdata    = testdata';

[freal, ~]  = prob.evaluate_l([], testdata);
plot(testdata, freal, 'b', 'LineWidth', 1);hold on;

% (2) plot archive
scatter(arc_xl, arc_fl, 80, 'ko');

% (3) plot local region
scatter(local_xl, local_fl, 90, 'rx');

% (4) check krging rebuilt
local_ub    = max(local_xl, [], 1);
local_lb    = min(local_xl, [], 1);

local_test  = linspace(local_lb, local_ub, 100);
local_test  = local_test';
[fpred, ~]  = surrogate_predict(local_test, krg, arc_obj);

% (5) denormalize prediction
fpred       = denormzscore(local_fl, fpred);

plot(local_test, fpred, 'r', 'LineWidth', 2);

% (6) plot local search result
scatter(localnewx,localnewf, 120, 'b*');

pause(1);
end

function plotlocalKprocess2D(fighn, prob, local_xl, local_fl, krg, krgc, ...
    localnewx, localnewf, arc_obj, arc_con, xu, varargin)
if length(krgc) > 1
    fprintf('not for multiple constraints \n');
    return;
end


clf(fighn);
set(fighn,'Position',[1000 10 1000 800])
% (1) create 2d test data and true local function
nt                  = 100;
local_bl            = min(local_xl, [], 1); % prob.xl_bl;
local_bu            = max(local_xl, [], 1); % prob.xl_bu;

if sum(local_bl >= prob.xl_bl) == size(local_xl, 2) % local plot
    delta_bl            = local_bl - prob.xl_bl;
    delta_bu            = prob.xl_bu - local_bu;
    x1_tst              = linspace(local_bl(1)-  0.5 * delta_bl(1), local_bu(1) + 0.5 * delta_bu(1), nt);
    x2_tst              = linspace(local_bl(2)- 0.5 * delta_bl(2), local_bu(2) + 0.5 * delta_bu(2), nt);
else % main process plot
    x1_tst              = linspace(local_bl(1), local_bu(1), nt);
    x2_tst              = linspace(local_bl(2), local_bu(2), nt);
    
end

[msx1, msx2]        = meshgrid(x1_tst, x2_tst);
f                   = zeros(nt, nt);
c1                  = zeros(nt, nt);

for i =1 : nt
    for j = 1: nt
        if ~isempty(krgc)
            [f(i, j),c1(i, j)] ...
                = prob.evaluate_l(xu, [msx1(i, j), msx2(i, j)]);
        else
            f(i, j) ...
                = prob.evaluate_l(xu, [msx1(i, j), msx2(i, j)]);
        end
    end
end
subplot(2, 2, 1)
% contour(msx1, msx2, f); hold on;
surf(msx1, msx2, f, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;
title('real function');
colorbar
scatter3(localnewx(1), localnewx(2),  localnewf, 80, 'r', 'filled' ); hold on;
scatter3(local_xl(:, 1), local_xl(:, 2), local_fl,  60, 'b');
xlabel('xl-1'); ylabel('x1-2');
xlim([prob.xl_bl(1), prob.xl_bu(1)]);
ylim([prob.xl_bl(2), prob.xl_bu(2)]);


subplot(2, 2, 2)
contour(msx1, msx2, c1); hold on;
title('real constraints');
colorbar
scatter(localnewx(1),localnewx(2),  80, 'r', 'filled' );
scatter(local_xl(:, 1), local_xl(:, 2),  80, 'kx' );
xlabel('xl-1'); ylabel('x1-2');
xlim([prob.xl_bl(1), prob.xl_bu(1)]);
ylim([prob.xl_bl(2), prob.xl_bu(2)]);

%---------------------------------------------------------------------
% create prediction map
local_bl            = min(local_xl,[], 1);
local_bu            = max(local_xl,[], 1);

% nt                  = 200;
fp                  = zeros(nt, nt);
cp                  = zeros(nt, nt);
x1_tst              = linspace(local_bl(1), local_bu(1), nt);
x2_tst              = linspace(local_bl(2), local_bu(2), nt);
[msx1, msx2]        = meshgrid(x1_tst, x2_tst);
for i =1 : nt
    for j = 1: nt
        [fp(i, j), ~] = surrogate_predict( [msx1(i, j), msx2(i, j)], krg, arc_obj);
        fp(i, j)      = denormzscore(local_fl, fp(i, j));
        if ~isempty(krgc)
            [cp(i, j), ~]...
                = surrogate_predict( [msx1(i, j), msx2(i, j)], krgc, arc_con);
        else
            
            if ~isempty (varargin) % local search
                cp = fp;
                
            else % main process
                ynorm_min         = min(arc_obj.muf);
                cp(i, j)          = EIM_evaldaceUpdate([msx1(i, j), msx2(i, j)], ynorm_min,  arc_obj, krg);
                cp(i, j)          = -cp(i, j);
            end
            
        end
        
    end
end
trgidx      = my_ismember(arc_obj.x, local_xl);
a           = 1: size(local_xl, 1);
a(trgidx)   = [];
tstx        = local_xl(a, :);
tstf        = local_fl(a, :);


subplot(2, 2, 3);
% contour(msx1, msx2, fp); hold on;


surf(msx1, msx2, fp, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;
colorbar
title('surrogate function');
xlabel('xl-1'); ylabel('x1-2');
scatter3(localnewx(1),  localnewx(2), localnewf,  80, 'r', 'filled' );
arcmuf = denormzscore(local_fl, arc_obj.muf);
scatter3(arc_obj.x(:, 1), arc_obj.x(:, 2),arcmuf,  60, 'b', 'filled');
scatter3(tstx(:, 1), tstx(:, 2),tstf,  60, 'b');


subplot(2, 2, 4);

if isempty(krgc)
    

    
    if isempty (varargin) % main process
        surf(msx1, msx2, cp, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;
    
        colorbar
        ynorm_min       = min(arc_obj.muf);
        eif                     = EIM_evaldaceUpdate(localnewx, ynorm_min,  arc_obj, krg);
        eif                     = -eif;
        
        scatter3(localnewx(1),  localnewx(2), eif,  80, 'r', 'filled' ); hold on;
        title('EIM landscape');
    else
       
        surf(msx1, msx2, cp, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;
    
        colorbar
        title('local surrogate')
    end
else
    contour(msx1, msx2, cp); hold on;
    colorbar
    title('surrogate constraint1');
    scatter(localnewx(1),  localnewx(2), 80, 'r', 'filled' );
    scatter(local_xl(:, 1), local_xl(:, 2),  80,'kx' );
end
xlabel('xl-1'); ylabel('x1-2');
pause(0.5);

end

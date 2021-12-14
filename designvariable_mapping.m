% Bilevel mapping
% in this experiment I would like to test the possibility of using mapping of
% xu to xl as starting point for lower level landscape
% (1) Generate a population of 50
% (2) Run global search on lower level
% (3) Create the bilevel mapping
% (4) Predict on new landscape
% (5) Show plot


clear;
close all;



problem_folder = strcat(pwd, '/TP3');
addpath(problem_folder);
problem_folder = strcat(pwd, '/globalsolver');
addpath(problem_folder);
problem_folder = strcat(pwd, '/ND_Sort');
addpath(problem_folder);
problem_folder = strcat(pwd, '/OODACE');
addpath(problem_folder);

prob = smd4mp(1, 1, 1);

rng(1, 'twister');
% (1) 
nu = 50;
xu = lhsdesign(nu, prob.n_lvar, 'criterion', 'maximin', 'iterations', 1000);
xu = repmat(prob.xl_bl, nu, 1) + repmat((prob.xl_bu - prob.xl_bl), nu, 1) .* xu;

% (2)
xl = [];
param.gen = 9;
param.popsize = 50;


for i = 1:nu
    funh_obj = @(x)lower_objective(prob, xu(i, :), x);
    funh_con = @(x)lower_constraint(prob, xu(i, :), x);
    
    [best_x, best_f, best_c, archive ,~] = gsolver(funh_obj, prob.n_lvar, prob.xl_bl, prob.xl_bu, [], funh_con, param, 'visualize', false);
    xl = [xl; best_x];
end

% (3)
mdls = {};
for i = 1 : prob.n_lvar
    mdls{i} = oodacefit(xu, xl(:, i));
end

ns = 7;
xu_tst = lhsdesign(ns, prob.n_lvar, 'criterion', 'maximin', 'iterations', 1000);
xu_tst = repmat(prob.xl_bl, ns, 1) + repmat((prob.xl_bu - prob.xl_bl), ns, 1) .* xu_tst;


% (4) check prediction
k = 7;
xltst = [];
for i = 1:prob.n_lvar
    tmp = mdls{i}.predict(xu_tst(k, :));
    xltst = [xltst, tmp]; 
end

% (5) plot neighbour
dist = pdist2(xu_tst(k, :) , xu);                              % this xu is upper level new infill xu, not added into archive_xu
[~, idx] = sort(dist);
xu_neighbour = xu(idx(1), :);
xl_neighbour = xl(idx(1), :);

% (5) plot
plotlower(prob, prob.xl_bl, prob.xl_bu, xu_tst(k, :), xltst, xu_neighbour, xl_neighbour);

function plotlower(prob, lb, ub, xu, xl_pred, xu_neighbour, xl_neighbour)
fighn = figure(1);

nt = 100;
x1_tst              = linspace(lb(1), ub(1), nt);
x2_tst              = linspace(lb(2), ub(2), nt);
[msx1, msx2]        = meshgrid(x1_tst, x2_tst);
msx11 = msx1(:);
msx22 = msx2(:);
xl  = [msx11, msx22];
xum = repmat(xu, nt*nt, 1);
fl  = prob.evaluate_l(xum, xl);
fl  = reshape(fl, [nt, nt]);

subplot(1, 2, 1);
surf(msx1, msx2, fl,'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;
fl_pred = prob.evaluate_l(xu, xl_pred);
fl_neighbour = prob.evaluate_l(xu, xl_neighbour);


scatter3(xl_pred(1), xl_pred(2), fl_pred,80, 'red', 'filled'); hold on;
scatter3(xl_neighbour(1), xl_neighbour(2), fl_neighbour, 40, 'blue', 'filled'); hold on;


subplot(1, 2, 2);
xum = repmat(xu_neighbour, nt*nt, 1);
fl  = prob.evaluate_l(xum, xl);
fl  = reshape(fl, [nt, nt]);
surf(msx1, msx2, fl,'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;
fl_neighbour = prob.evaluate_l(xu_neighbour, xl_neighbour);
scatter3(xl_neighbour(1), xl_neighbour(2), fl_neighbour, 80, 'blue', 'filled'); hold on;


pause(1);

end






function obj = lower_objective(prob, xu, xl)
obj = prob.evaluate_l(xu, xl);
end

function con = lower_constraint(prob,xu, xl)
con = [];
end






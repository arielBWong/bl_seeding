function [bestx, bestf, bestc, archive, external_return] = ego_solver(funh_obj, num_xvar, lb, ub, initmatrix, funh_con, param, varargin)
% this function is a wrapper to EGO optimization
% processing single objective only
% input
%       funh_obj : function handle to objective function
%       num_xvar : number design variables
%       lb: upper bound of design varibles
%                1d row array
%       up: lower bound of design variables
%       initmatrix:  partial population to be embeded in
%                           initial population
%       funh_con : function handle to constraint functions
%       param : structure specifying kriging parameters(param.maxFE; param.initsize)
%       varargin : don't know what to do with it yet
% output:
%       bestx : global search results of design variables
%       bestf : global search results of objective values
%       bestc : global search results of constraints
%                       for constraint problems, if no feasible is found, return least infeasible one
%
%--------------------------------------------------------------------------
% create argument parser
p = inputParser;
addRequired(p, 'funh_obj');
addRequired(p, 'num_xvar');
addRequired(p, 'lb');
addRequired(p, 'ub');
addRequired(p, 'initmatrix');
addRequired(p, 'funh_con');
addRequired(p, 'param');
addParameter(p, 'externalfunction', []);
addParameter(p, 'visualize', false);
addParameter(p, 'ea_param', []);
parse(p, funh_obj, num_xvar, lb, ub, initmatrix, funh_con, param, varargin{:});

%-------------------
funh_obj = p.Results.funh_obj;
num_xvar = p.Results.num_xvar;
lb = p.Results.lb;
ub = p.Results.ub;
initmatrix = p.Results.initmatrix;
funh_con = p.Results.funh_con;
param = p.Results.param;
external_funh = p.Results.externalfunction;
visualize = p.Results.visualize;
ea_param = p.Results.ea_param;
%-------------------

external_return = [];
bestx = NaN;
bestf = NaN;
bestc = NaN;

%--------start--------------------
n_init = size(initmatrix, 1);
n_rest = param.initsize  - n_init;

trgx = repmat(lb, n_rest, 1) + repmat(ub - lb, n_rest, 1) .* lhsdesign(n_rest, num_xvar);
trgx = [trgx; initmatrix];
trgf =  funh_obj(trgx);
if visualize
    fighn = figure(1);
end

%--- train model
mdl = oodacefit(trgx, trgf);

% infill process
while size(trgx, 1) < param.maxFE
    % searching for next infill using believer
    
    ea_param.gen = 100;                    % design problem of gsolver, it should be 20
    ea_param.popsize = 100;
    
    infill_obj = @(x)infill_objective(mdl, x);
    infill_con = @(x)infill_constraints(x);
    
    % initmatrix is set early in this method;
    [next_x, ~, ~, ~, ~] = gsolver(infill_obj, num_xvar, lb, ub, [], infill_con, ea_param, 'visualize', false);
    
    % add to trg data
    next_f = funh_obj(next_x);
    
    % place visualization before update  trg
    if visualize
        plot2d_infill(fighn, lb, ub, funh_obj, mdl, trgx, trgf, next_x);
    end
    
    % training data extension needs distance check
    [trgx, trgf] = trgdata_extension_withDistCheck(trgx, trgf, next_x, next_f, lb, ub);
    
    %-------------
    mdl = oodacefit(trgx, trgf);
end

if visualize
   close(fighn);
end

% select best solution
[~, idx] = sort(trgf);
bestx = trgx(idx(1), :);
bestf = trgf(idx(1), :);

% only deal with single objective problems
bestc = [];
n = 1:size(trgx, 1);
archive.sols = [n', trgx, trgf];



end

function f = infill_objective(mdl, x)
[f, ~] = mdl.predict(x);
end

function c = infill_constraints(x)
c = [];
end

function plot2d_infill(fighn, lb, ub, funh_obj, mdl, trgx, trgf, infillx)

clf(fighn);

% (1) generate test process
nt = 100;
tst1 = linspace(lb(1), ub(1), nt);
tst2 = linspace(lb(2), ub(2),  nt);
[msx1, msx2] = meshgrid(tst1, tst2);
msx11 = msx1(:);
msx22 = msx2(:);

truef = funh_obj([msx11, msx22]);
truef = reshape(truef, [nt, nt]);

subplot(1, 2, 1);
surf(msx1, msx2, truef, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;

infillf = funh_obj(infillx);
scatter3(infillx(1), infillx(2), infillf,  80, 'r', 'filled' );
scatter3(trgx(:, 1), trgx(:, 2), trgf, 80, 'g', 'filled' );


subplot(1, 2, 2);
predf = mdl.predict([msx11, msx22]);
predf = reshape(predf, [nt, nt]);

surf(msx1, msx2, predf, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
scatter3(trgx(:, 1), trgx(:, 2), trgf, 80, 'g', 'filled' );
infillf = mdl.predict(infillx);
scatter3(infillx(1), infillx(2), infillf,  80, 'r', 'filled' );

pause(0.5);
end

function [trgx, trgf] = trgdata_extension_withDistCheck(trgx, trgf, newx, newf, lb, ub)
trgx_norm = (trgx - lb) ./ (ub - lb);
newx_norm = (newx - lb) ./ (ub - lb);

[dist, I] = pdist2(trgx_norm, newx_norm, 'euclidean', 'Smallest', 1);

if dist < 1e-6 % new point is too close
    if  newf < trgf(I)  % new point is better
        trgx(I, :) = [];
        trgf(I, :) = [];
        trgx = [trgx; newx];
        trgf = [trgf; newf];
    end
    % else no change
else
    trgx = [trgx; newx];
    trgf = [trgf; newf];
end




end
function [bestx, bestf, bestc, archive, external_return] = gsolver(funh_obj, num_xvar, lb, ub, initmatrix, funh_con, param, varargin)
% This function is a wrapper on methods in global optimization/minimization
% folder. Main process is nsga2, child population is generated by a DE operator
% Be aware, this method only handle inequality constraints
% Don't use on equality constraints
% input:
%       funh_obj : function handle to objective function
%       num_xvar : number design variables
%       lb: upper bound of design varibles
%                1d row array
%       up: lower bound of design variables
%       initmatrix:  partial population to be embeded in
%                           initial population
%       funh_con : function handle to constraint functions
%       param : structure specifying ea parameters(param. popsize; param.gen)
%       varargin : additional variables for dealing with bilevel problems
% output:
%       bestx : global search results of design variables   (best value or nd front)
%       bestf : global search results of objective values   (best value or nd front)
%       bestc : global search results of constraints
%                       for constraint problems, if no feasible is found, return least infeasible one
% 
%--------------------------------------------------------------------------
% create argument parser

p = inputParser;
addRequired(p,   'funh_obj');
addRequired(p,   'num_xvar');
addRequired(p,   'lb');
addRequired(p,   'ub');
addRequired(p,   'initmatrix');
addRequired(p,   'funh_con');
addRequired(p,   'param');
addParameter(p,  'externalfunction', []);
addParameter(p,  'visualize', false);
addParameter(p, 'output_selection', []);
addParameter(p, 'extra_infoprob', []);
parse(p, funh_obj, num_xvar, lb, ub, initmatrix, funh_con, param, varargin{:});
%-------

funh_obj= p.Results.funh_obj;
num_xvar = p.Results.num_xvar;
lb = p.Results.lb;
ub = p.Results.ub;
initmatrix= p.Results.initmatrix;
funh_con = p.Results.funh_con;
param = p.Results.param;
external_funh= p.Results.externalfunction;
visualize = p.Results.visualize;
output_selection = p.Results.output_selection;
extra_infoprob = p.Results.extra_infoprob;

%-----------
external_return = [];
bestx = NaN;
bestf = NaN;
bestc = NaN;

if visualize
    f1 = figure(2);
    obj = VideoWriter('moving.avi');
    obj.Quality= 100;
    obj.FrameRate = 25;
    open(obj);
end


% Initialization
[pop,archive] = initialize_pop(funh_obj, funh_con, num_xvar, lb, ub, initmatrix, param);

if ~isempty(external_funh)
    external_funh(pop);
end

%%%%%%%%%%%%%%%%%%%%%
if visualize && size(lb, 2) == 1
     plot1d(f1, lb, ub, param)
end

if visualize && size(lb, 2)== 2 && isempty(external_funh) 
     plot2d(f1, lb, ub, param, pop, funh_obj);
end

if visualize && size(lb, 2) == 2 && ~isempty(external_funh)
     plot2dupper(f1, lb, ub, pop);
end
%%%%%%%%%%%%%%%%%%%%


gen=1;
while gen <= param.gen
    % Recombination
    % child.X=generate_child(lb, ub, pop, param);
    child.X = generate_child2(lb, ub, pop, param);
    
    % Evaluate and Order(need to test varargin)
    % [pop, archive]= evaluate_order(pop, archive, funh_obj, funh_con, child.X, gen, param);
    
     [pop, archive]= evaluate_order2(pop, archive, funh_obj, funh_con, child.X, gen, param, 'prob', extra_infoprob);
    
    % Reduce 2N to N
    [pop]=reduce_pop(pop,param.popsize);
    
    % for external function
    if ~isempty(external_funh)
        external_funh(pop);
    end
  
    %%%%%%%%%%%%%%%%%%%%%%%%%
    if visualize && size(lb, 2) == 1
       plot1d(f1, lb, ub, param); 
    end
    
    if  visualize && size(lb, 2) == 2&& isempty(external_funh) 
        plot2d(f1, lb, ub, param, pop, funh_obj);
%         frame = getframe(gcf);
%          for i = 1: 20
%              writeVideo(obj, frame);
%          end
        
    end
    
     if visualize && size(pop.F, 2) == 3
        plotMO(f1, pop.F, varargin{2});
     end
    
     if visualize && size(lb, 2) == 2 && ~isempty(external_funh)
         plot2dupper(f1, lb, ub, pop);
%          frame = getframe(gcf);
%          for i = 1: 20
%              writeVideo(obj, frame);
%          end
         
     end

    %%%%%%%%%%%%%%%%%%%%%%%%
    gen = gen+1;
    % disp(gen);
end



if visualize
    close(f1); 
    obj.close();
end


% use archive to save last pop_x
archive.pop_last =pop;

num_obj = size(pop.F, 2);
num_con = size(pop.C, 2);

if num_obj == 1
    bestx = pop.X(1, :);
    bestf = pop.F(1, :);
    if ~isempty(pop.C)
        bestc = pop.C(1, :);
    else
        bestc = [];
    end
    return
end

% return feasible nd front
if num_obj > 1
    if ~isempty(pop.C)                                           % constraint problem
        pop.C(pop.C<=0) = 0;
        fy_ind = sum(pop.C, 2) ==0;
    else
        fy_ind =true(param.popsize ,1);                    % unconstraint problem
    end
end

% feasible exists for mo problem
if sum(fy_ind) > 0
    [fronts, ~, ~] = nd_sort(pop.F, find(fy_ind));
    if isempty(output_selection)
        bestf = pop.F(fronts(1).f, :);
        bestx = pop.X(fronts(1).f, :);
        
        
        if ~isempty(pop.C)
            bestc = pop.C(fronts(1).f, :);
        else
            bestc = [];
        end
    else
        [bestx, bestf, bestc ] = output_selection(fronts, pop);
        
    end
else
    % no feasible solution in final population
    % return least infeasible on
    sumc= sum(pop.C, 2);
    [~, id] = sort(sumc);
    bestf = pop.F(id(1),  :);
    bestx = pop.X(id(1),  :);
    bestc = pop.C(id(1),  :);
end
end


function plotMO(fighn, pop, pf )
clf(fighn);
plot3(pf(:, 1), pf(:, 2), pf(:, 3), 'go'); hold on;
plot3(pop(:, 1), pop(:, 2), pop(:, 3), 'ro');
grid on;
pause(0.1);
end

function plot1d(fighn, lb, ub, param, pop)
    clf(fighn);
    xlim([lb, ub]);
    testdata        = linspace(lb, ub, 10000);
    testdata        = testdata';
    
    ynorm_min       = param.ynorm_min;
    arc_obj         = param.arc_obj ;
    krg             = param.krg  ;
    fit             = EIM_evaldaceUpdate(testdata, ynorm_min,  arc_obj, krg);
    fit             = -fit;
    plot(testdata, fit, 'g--'); hold on;
    scatter(pop.X, -pop.F, 50, 'filled', 'r');
    xlim([lb, ub]);
    drawnow;
    
end

function plot2d(fighn, lb, ub, param, pop, funh_obj)
clf(fighn);
nt                  = 100;

cp                  = zeros(nt, nt);
x1_tst              = linspace(lb(1), ub(1), nt);
x2_tst              = linspace(lb(2), ub(2), nt);
[msx1, msx2]        = meshgrid(x1_tst, x2_tst);



for i =1 : nt
    for j = 1: nt
        cp(i, j)    = funh_obj([msx1(i, j), msx2(i, j)]);
        % cp(i, j)    = -cp(i, j);
    end
end
surf(msx1, msx2, cp, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;
scatter3(pop.X(:, 1), pop.X(:, 2), pop.F,  80, 'r', 'filled' ); hold on;

pause(0.5);
drawnow;

end

function plot2dupper(fighn, lb, ub, pop)
global prob
clf(fighn);
nt                  = 100;

cp                  = zeros(nt, nt);
x1_tst              = linspace(lb(1), ub(1), nt);
x2_tst              = linspace(lb(2), ub(2), nt);
[msx1, msx2]        = meshgrid(x1_tst, x2_tst);
msx11 = msx1(:);
msx22 = msx2(:);
xu = [msx11, msx22];
xl = prob.get_xlprime(xu);

fu = prob.evaluate_u(xu, xl);
fu = reshape(fu, [nt, nt]);

surf(msx1, msx2, fu, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;
scatter3(pop.X(:, 1), pop.X(:, 2), pop.F,  80, 'r', 'filled' ); hold on;


pause(0.5);

end

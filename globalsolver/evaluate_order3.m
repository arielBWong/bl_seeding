function [pop, archive, termination_flag]= evaluate_order3(pop, archive, funh_obj, funh_con, cx, gen, param, varargin)
% This is a customized sorting scheme
% There are two sorting steps. The first sorting is to check whether
% children solution is based on
% input
%           pop : population of previous generation
%           archive : record of all evolutionary process
%           funh_obj : function handle of objective function
%           funh_con: function handle of constraints
%           cx: child population generated from last step
%           gen: current generation
%           param: evolution parameter (gen popsize)
% output
%           pop: extended and sorted population
%           archive: extended archive
%-------------------------
p = inputParser;
addRequired(p, 'pop');
addRequired(p, 'archive');
addRequired(p, 'funh_obj');
addRequired(p, 'funh_con');
addRequired(p, 'cx');
addRequired(p, 'gen');
addRequired(p, 'param');
addParameter(p, 'prob', []);
parse(p, pop, archive, funh_obj, funh_con, cx, gen, param, varargin{:});

%--------------------
pop = p.Results.pop;
archive = p.Results.archive;
funh_obj = p.Results.funh_obj;
funh_con = p.Results.funh_con;
cx = p.Results.cx;
gen = p.Results.gen;
param = p.Results.param;
prob = p.Results.prob;
%--------------------
child.X = cx;

% This is NSGA-II
output = funh_obj(child.X);
termination_flag = false;

if isstruct(output)     % cumstomized
    child.F = output.f;
    child.A = output.addon;
    child.Mdl = output.mdl;
    child.trgdata = output.trgdata;
    child.switch_lls = output.lower_searchSwitchFlags;
    child.LLcount = output.LLcount;

    if output.termination_flag == true
        termination_flag = true;
        child.X = child.X(1:size(child.A, 1), :);
    end
else
    child.F = output;
    child.A = [];
    child.Mdl = {};
    child.trgdata = {};
    child.switch_lls = [];
    child.LLcount = [];
end


child.C = funh_con(child.X);
numcon = size(child.C, 2);
archive.sols = [archive.sols;[repmat(gen, size(child.X, 1), 1), child.X, child.F, child.C, child.A]];

% Appending X F C to pop
pop.X = [pop.X; child.X];
pop.F = [pop.F; child.F];
pop.C = [pop.C; child.C];
pop.A = [pop.A; child.A];   % add on
pop.Mdl = [pop.Mdl, child.Mdl];
pop.trgdata = [pop.trgdata, child.trgdata];
pop.switch_lls = [pop.switch_lls; child.switch_lls];
pop.LLcount = [pop.LLcount; child.LLcount];

% [pop.F, pop.X, pop.C, pop.A, pop.Mdl, pop.trgdata] = pop_sort(pop.F, pop.X, pop.C, pop.A, pop.Mdl, pop.trgdata);
pop = pop_sortCompact(pop);

% check whether the first solution is local search
% means local search
if ~isempty(pop.switch_lls)

    if pop.switch_lls(1) == 1        
        % recalcualte the  lower level re-evaluation
        xu = pop.X(1, :);
        initmatrix = pop.trgdata{1};
        initmatrix = initmatrix(:, 1:end-1);
        num_xvar = size(pop.A, 2);
        funh_objego = @(x)objective_func(prob, xu, x);
        funh_conego = @(x)constraint_func(prob, xu, x);
        select_hn = @output_selection;
        param.maxFE = 150;
        param.initsize = 50;
        [best_x, best_f, best_c, archive_search] = ego_solver(funh_objego, num_xvar, prob.xl_bl, prob.xl_bu, initmatrix, funh_conego, param, 'visualize', false,...
            'infill', 3, 'gsolver_outputselect', select_hn, 'initmatrix_asExtra', true);        
        
        % follow local search
        local_FE = 50;
        % attach a local search
        [xsol, ~, history, output] = lowerlevel_fmincon(best_x, local_FE, prob.xl_bl, prob.xl_bu,funh_objego, funh_conego);
        localdata = [history.x, history.fval];
        
        trgdata = [archive_search.sols(:, 2: prob.n_lvar + 1), archive_search.sols(:, prob.n_lvar + 2: end)];
        trgdata = [trgdata; localdata];
              
        pop.switch_lls(1) = 0;
        pop.F(1) = prob.evaluate_u(xu, xsol);
        pop.C = [];           % not considering constraints
        pop.A(1, :) =  xsol;
        pop.trgdata(1) = {trgdata};   
        pop.LLcount(1) = pop.LLcount(1) + param.maxFE + output.funcCount;

        global lower_eval
        lower_eval = lower_eval + param.maxFE + output.funcCount;


        % Can I change it to more efficient in terms of convergence?
        % resort
        pop =  pop_sortReval(pop);
        
    end
end

% % dealing with feasibility on lower level
% if ~isempty(varargin)
%     fprintf('in evaluate_order, customized, but test needed\n');
%     prob = varargin{1};
%     xu_g =  varargin{2};
%     xl_g = varargin{3};
%     pop.F = pop_llfeasicheck(pop.X, pop.F, xu_g, xl_g, prob);
% end

end

function [pop] = pop_sortReval(pop)
%

num_pop = size(pop.F, 1);

% This  function
indx = 1: num_pop;
indx_global  = pop.switch_lls==0;
indx_local = pop.switch_lls == 1;

indx_globalnum = indx(indx_global);
indx_localnum = indx(indx_local);

popF_global = pop.F(indx_global, :);
popF_local = pop.F(indx_local, :);

%
[~, idx_global]  = sort(popF_global);
[~, idx_local] = sort(popF_local);

ids= [indx_globalnum(idx_global), indx_localnum(idx_local)];

pop.F = pop.F(ids, :);
pop.X = pop.X(ids, :);

if ~isempty(pop.A)
    pop.A =  pop.A(ids,:);
end

if ~isempty(pop.Mdl)
    pop.Mdl = pop.Mdl(ids);   % cell array
end

if ~isempty(pop.trgdata)
    pop.trgdata = pop.trgdata(ids);   % cell array
end

if ~isempty(pop.switch_lls)
    pop.switch_lls = pop.switch_lls(ids);   % cell array
end

if ~isempty(pop.LLcount)
    pop.LLcount = pop.LLcount(ids, :);
end


end


function f = objective_func(prob, xu, xl)
n = size(xl, 1);
xu = repmat(xu, n, 1);

f = prob.evaluate_l(xu, xl);
end

function [c, ce] = constraint_func(prob, xu, xl)
c = [];
ce = [];
end

function[mf] = pop_llfeasicheck(xu, fu, xu_g, xl_g, probh)
% this function considers bilevel situation, where lower level is infeasible
%  it collects all upper level functions that get evaluated (it is infact all past population)
% uses its highest values (minimization problem)  to punish ll infeasible
% solutions
% usage
%   input
%  xu : current population
%  fu:  current population objectives to be re-accessed
%  xu_g : global variable passed from outside
%  xl_g :  global variable passed from outside
% probh: problem instance
%   output
%   : modified current population on fu
% * test required
%-----------------------------------------------------------------

% get feasibility from global saves
% since global xu xl are already evaluated
% here calculating fc is not counted for algorithm clearity in
% compatible with addons for gslover, sacrifices efficiency for
% compatibility
workdir = pwd;
idcs = strfind(workdir, '\');
upperfolder = workdir(1: idcs(end)-1);
addpath(upperfolder);   % llfeasi_modify exists in upper folder

num_pop = size(fu, 1);
xl = [];
for i = 1:num_pop
    match =  check_exist(xu(i, :), xu_g, xl_g);
    if isempty(xl)
        error('logically wrong, by far, every xu in population should have xl in global variable');
    end
    xl = [xl; match];
end

[~, fc] = probh.evaluate_l(xu, xl);
if ~isempty(fc)
    % in order to reuse function llfeasi_modify
    % construct arguments for llfeasi_modify
    [~, fc_g]  = probh.evaluate_l(xu_g, xl_g);
    [fu_g, ~] = probh.evaluate_u(xu_g, xl_g);
    num_con = size(fc_g, 2);
    num_g = size(fc_g, 1);
    % (1) feasi_flag of all solutions
    fc_all = [fc_g; fc];
    fc_all(fc_all<=0) = 0;
    feasi_flag = sum(fc_all, 2)==0;
    % (2) contruct all current evaluated fus
    fu_all = [fu_g; fu];
    % (3) process fu one by one use index in all solutions
    for ii = num_g + 1: num_g + num_pop
        fu_all = llfeasi_modify(fu_all, feasi_flag, ii);
    end
    mf = fu_all(num_g + 1: num_g + num_pop, :);
    
else
    mf = fu;
end
%----------------------------------------------
end


function xl = check_exist(xu, xu_g, xl_g)
%  this function should can be extracted
% as it is also used in ul_llea m file
%--------------------------------------
xu = round(xu, 10);
if isempty(xu_g)
    xl = [];
else
    xu_g = round(xu_g, 10);
    diff = xu_g - xu;
    ind = sum(diff, 2) == 0;
    
    if sum(ind)>0           % should be only 1, if exists
        % fprintf('found xu in global save %d\n', sum(ind));
        xl = xl_g(ind, :);
        if sum(ind)>1
            error('there cannot be repeat x more than once')
        end
    else
        xl = [];
    end
end
end


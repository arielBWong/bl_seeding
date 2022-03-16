function[match_xl, mdl, trgdata, lower_searchSwitchFlag] = llmatch_egoEvaluation(upper_xu, prob, xl_probe, varargin)
% method of searching for a match xl for xu.
% Problem(Prob) definition require certain formation for bilevel problems
% evaluation method should be of  form 'evaluation_l(xu, xl)'
%--------------------------------------------------------------------------
% distribute parameters
% parse input
%-------
p = inputParser;
addRequired(p, 'upper_xu');
addRequired(p, 'prob');
addParameter(p, 'visualization', false);
addParameter(p, 'seeding_only', false);      % conduct seeding controller
addParameter(p, 'archive_xu', []);               % upper level all evaluated xu
addParameter(p, 'archive_xl', []);                % lower level matching xl, might larger than xu due to asynsic update
addParameter(p, 'seeding_strategy', 1);     % 1 nearest neighbour 2 cokrging starting
addParameter(p, 'threshold', 0.9);              % 1 nearest neighbour 2 cokrging starting
parse(p, upper_xu, prob, xl_probe, varargin{:});
% -----------------------------------

xu = p.Results.upper_xu;
visualization = p.Results.visualization;
seeding_only = p.Results.seeding_only;
archive_xu = p.Results.archive_xu;
archive_xl = p.Results.archive_xl;
prob = p.Results.prob;
seeding_strategy = p.Results.seeding_strategy;
thr = p.Results.threshold;
%------------------------------------------------

lower_searchSwitchFlag = 0;
global lower_eval
global lower_trg
%-------------------------------------------------

cokrg_samplesize = 20;
mdl = 'position_holder';
trgdata = 'position_holder';
initmatrix = []; %

funh_obj = @(x)objective_func(prob, xu, x);
funh_con = @(x)constraint_func(prob, xu, x);

param.maxFE = 150;
param.initsize = 50;

%-- lower level initialization
% include seeding xl
if ~isempty(archive_xu) && seeding_only % first generation on the upper level use global search
    % ideal_xl = prob.get_xlprime(xu);
    maxFE = 200;    % maxFE for local search

    if seeding_strategy == 1  % closest landscape
        [~, idx]= pdist2(archive_xu, xu, 'euclidean', 'Smallest', 1);  % this xu is upper level new infill xu, not added into archive_xu
        close_optxl =  archive_xl(idx, :);
        [match_xl, ~, history, output] = lowerlevel_fmincon(close_optxl, maxFE, prob.xl_bl, prob.xl_bu, funh_obj, funh_con);

        lower_eval = lower_eval + output.funcCount;
        trgdata = [history.x, history.fval];
        lower_searchSwitchFlag = 1;
        return;
    end

    if seeding_strategy == 2    % test neighbour
        % this method use cokriging sample to determine a starting point
        % and local search
        % cokriging sample both consider in std and outside std
        [expensive_x, expensive_f, cheap_x, cheap_f, correlation] = ...
            cokrg_trainingSampleExtension(xu, prob, archive_xu, archive_xl,  lower_trg, cokrg_samplesize);

        [neighbour_optxu, neighbour_optxl, ~] = retrieve_neighbour(xu, lower_trg, archive_xu, archive_xl);
        neighbour_optfl = prob.evaluate_l(neighbour_optxu, neighbour_optxl);
        current_fl = prob.evaluate_l(xu, neighbour_optxl);

        cokrg_trg.expensive_x = [expensive_x; neighbour_optxl];
        cokrg_trg.expensive_f = [expensive_f; current_fl];
        cokrg_trg.cheap_x = [cheap_x; neighbour_optxl];
        cokrg_trg.cheap_f = [cheap_f; neighbour_optfl];

        if correlation > thr
            % expensive x went through distance check
            maxFE = maxFE - size(cokrg_trg.expensive_x, 1);

            % expensive_x is counted inside cokrg_localsearch
            [match_xl, lower_eval, trgdata, localsearch_fail, cokrg_mdl] = cokrg_localsearch(xu, prob, cokrg_trg, lower_eval, maxFE);

            % if cokrging propose a point on boundary
            if localsearch_fail
                initmatrix = [cokrg_trg.expensive_x; match_xl];
                param.maxFE = param.maxFE - size(cokrg_trg.expensive_x, 1) - 2; % fail means two more evaluations in cokrg

                %  algorithm will continue to ea search part
            else
                lower_searchSwitchFlag = 1;

%                 % create test
%                 xl_prime = prob.get_xlprime(xu);
%                 fl_prime = prob.evaluate_l(xu, xl_prime);
%                 fl_pred = prob.evaluate_l(xu, match_xl);
%                 if abs(fl_prime - fl_pred) > 5
%                     fprintf('transfer  search fails \n');
%                     % visualize_globalSearch(prob, trgdata, xu, match_xl);
%                 end

                %----
                return;
            end
        else
            initmatrix = cokrg_trg.expensive_x;
            % initmatrix = [];
            param.maxFE = param.maxFE - size(cokrg_trg.expensive_x, 1);
        end
    end

    if seeding_strategy == 3    % test local search only
        % this method use cokriging sample to determine a starting point
        % and local search
        % cokriging sample both consider in std and outside std
        [expensive_x, expensive_f, ~, ~, correlation] = ...
            cokrg_trainingSampleExtension(xu, prob, archive_xu, archive_xl,  lower_trg, cokrg_samplesize);

        maxFE = maxFE - size(expensive_x, 1);
        if correlation > thr
            % retrieve ID
            [~,~, close_id] = retrieve_neighbour(xu, lower_trg, archive_xu, archive_xl);
            close_optxl =  archive_xl(close_id, :);
            close_optxu = archive_xu(close_id, :);
            [match_xl, ~, history, output] = lowerlevel_fmincon(close_optxl, maxFE, prob.xl_bl, prob.xl_bu, funh_obj, funh_con);

            % record training data and evaluation usage
            trgdata = [expensive_x, expensive_f;  ...
                history.x, history.fval];
            lower_eval = lower_eval + size(expensive_x, 1) + output.funcCount;
            lower_searchSwitchFlag = 1;


            if visualization
                neighbor_f = prob.evaluate_l(close_optxu, close_optxl);
                neighbor_optxl = prob.get_xlprime(close_optxu);
                neighbor_optf = prob.evaluate_l(close_optxu, neighbor_optxl);

                current_f = prob.evaluate_l(xu, match_xl);
                current_optxl = prob.get_xlprime(xu);
                current_optf = prob.evaluate_l(xu, current_optxl);

                %  if  abs(neighbor_f - neighbor_optf) < 0.5 &&  abs(current_f -current_optf) > 1
                %     plot2d_withNeighbor(xu, prob.xl_bl, prob.xl_bu, match_xl, close_optxl, prob, close_optxu)
                %  end
            end

            return;
        else
            initmatrix = expensive_x;
            % initmatrix = [];
            param.maxFE = param.maxFE - size(expensive_x, 1);
        end
    end
end


if visualization
    fighn            = figure(1);
end

% insert global search
% run to here means transfer fails or first generation

% param.gen = 19;                    % design problem of gsolver, it should be 20
% param.popsize = 50;
% lb = prob.xl_bl;
% ub = prob.xl_bu;
% num_xvar = prob.n_lvar;
% % initmatrix is set early in this method;
%
% [best_x, best_f, best_c, archive_search, ~] = gsolver(funh_obj, num_xvar, lb, ub, initmatrix, funh_con, param, 'visualize', false);

% apply believer Kriging

num_xvar = prob.n_lvar;

select_hn = @output_selection;
[best_x, best_f, best_c, archive_search] = ego_solver(funh_obj, num_xvar, prob.xl_bl, prob.xl_bu, initmatrix, funh_con, param, 'visualize', false,'infill', 3, 'gsolver_outputselect', select_hn);

% follow local search
local_FE = 50;
% attach a local search
[xsol, ~,history, output] = lowerlevel_fmincon(best_x, local_FE, prob.xl_bl, prob.xl_bu,funh_obj, funh_con);
localdata = [history.x, history.fval];

% recordings
lower_eval = lower_eval + param.maxFE - size(initmatrix, 1) + output.funcCount;  % initmatrix counted into  lower_eval in method cokrg_localsearch
match_xl = xsol;
trgdata = [archive_search.sols(:, 2: prob.n_lvar + 1), archive_search.sols(:, prob.n_lvar + 2: end)];
trgdata = [trgdata; localdata];


% % create test
% xl_prime = prob.get_xlprime(xu);
% fl_prime = prob.evaluate_l(xu, xl_prime);
% fl_pred = prob.evaluate_l(xu, xsol);
% if abs(fl_prime - fl_pred) > 5
%     fprintf('global search fails \n');
%     visualize_globalSearch(prob, trgdata, xu, xsol);
% end


% test visualization
%archive_xu = xu;
% archive_xl = match_xl;
% lower_trg{1} = trgdata;
% [expensive_x, expensive_f, cheap_x, cheap_f, correlation] = ...
%           cokrg_trainingSampleExtension(xu, prob, archive_xu, archive_xl,  lower_trg, cokrg_samplesize);


end



function f = denormzscore(trainy, fnorm)
[~, y_mean, y_std] = zscore(trainy);
f = fnorm * y_std + y_mean;
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

function visualize_globalSearch(prob, trgdata, xu, xsol)
fighn = figure(1);
nt = 100;
tst1 = linspace(prob.xl_bl(1), prob.xl_bu(1), nt);
tst2 = linspace(prob.xl_bl(2), prob.xl_bu(2), nt);

[tst11, tst22] = meshgrid(tst1, tst2);
tst11L = tst11(:);
tst22L = tst22(:);

xuu = repmat(xu, size(tst11L, 1), 1);

landscapefl = prob.evaluate_l(xuu, [tst11L, tst22L]);
landscapefl = reshape(landscapefl, [nt, nt]);

surf(tst11, tst22, landscapefl, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
scatter3(trgdata(:, 1), trgdata(:, 2), trgdata(:, 3),  80, 'g', 'filled' ); hold on;

fsol = prob.evaluate_l(xu, xsol);
scatter3(xsol(1), xsol(2), fsol, 100, 'r', 'filled' ); hold on;

close(fighn);


end


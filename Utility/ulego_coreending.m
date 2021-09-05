function ulego_coreending(xu, fu, fc, xl, prob, seed, n_up, n_low, method)
% fu is converted to one from archive
[xu_best, fu, cu, ~, index] = out_select(xu,  xl, prob);

name = prob.name;
fprintf('returned index %d for problem %s seed %d/n', index, name, seed);

xl_best = xl(index, :);
[fl, cl] = prob.evaluate_l(xu_best, xl_best);

cu(cu < 1e-6) = 0;
cl(cl < 1e-6) = 0;
% check feasibility
num_conu = size(cu, 2);
num_conl = size(cl, 2);
cu = sum(cu<=0, 2)==num_conu;
cl = sum(cl<=0, 2)==num_conl;
perf_record(prob, fu, cu, fl, cl, n_up, n_low, seed, method);
archive_record(xu, xl, prob, method, seed);

end

function archive_record(xu, xl, prob, method, seed)
num = prob.n_lvar;
savepath = strcat(pwd, '/result_folder/', prob.name, '_', num2str(num), '_', method);

n = exist(savepath);
if n ~= 7
    mkdir(savepath)
end

savename = strcat(savepath, '/xu_', num2str(seed),'.csv');
csvwrite(savename, xu);
savename = strcat(savepath, '/xl_', num2str(seed),'.csv');
csvwrite(savename, xl);

[fu, ~] = prob.evaluate_u(xu, xl);
[fl, ~] = prob.evaluate_l(xu, xl);

savename = strcat(savepath, '/fu_', num2str(seed),'.csv');
csvwrite(savename, fu);
savename = strcat(savepath, '/fl_', num2str(seed),'.csv');
csvwrite(savename, fl);

end


function [best_x,best_f, best_c, s, index] = out_select(xu,  xl, prob)
% this function select a starting point for a local search algorithm
% usage
% input:
%           x                       : 2d matrix of design variables,
%                                              (original range)
%           
% output:
%           best_x                  : one x instance with smallest objective value
%                                               or with smallest feasible objective value
%                                               or with smallest const value if no
%                                                               feasible solution exists
%           best_f                  : corresponding obj value of best_x
%           best_c                  : corresponding con value of best_x
%           s                       : flag indicate whether there is feasible solution
%           index                   : index of output x in archive
%--------------------------------------------------------------------------
[uf, uc] = prob.evaluate_u(xu, xl); % lazy fix
[lf, lc] = prob.evaluate_l(xu, xl); % lazy fix
c = [uc, lc];

if ~isempty(c)              		% constraint problem
     num_con = size(c, 2);
        index_c = sum(c <= 0, 2) == num_con;
        if sum(index_c) == 0        % no feasible, return f with smallest constraints
            sum_c = sum(c, 2);
            [~, i] = min(sum_c);
            best_x = xu(i, :);
            best_f = uf(i, :);
            best_c = uc(i, :);
            s = false;
            index = i;
        else 						% has feasible, return feasible smallest f
            feasi_ff = uf(index_c, :);
            feasi_x = xu(index_c, :);
            feasi_fc = uc(index_c, :);
            [~, i] = min(feasi_ff);
            best_x = feasi_x(i, :);
            best_f = feasi_ff(i, :);
            best_c = feasi_fc(i, :);
            s = true;
            [~, index] = ismember(best_x, xu, 'row');
            
        end
else                                % unconstraint problem
       [best_f, i] = min(uf);
        best_x = xu(i, :);
        s = true;
        best_c = [];
        index = i;
    
end


end
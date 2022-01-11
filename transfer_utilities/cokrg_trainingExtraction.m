function [expensive_x, cheap_x, cheap_f, correlation, close_optxu, close_optxl, global_xl, global_fl, sigma] = cokrg_trainingExtraction(xu , ...
    archive_xu, lower_trg, prob, global_halfflagsize, lower_xl, global_halfflag)
% this function identify close landscape for kriging (the one with global search)
% and identify training data for both Hx and Lx
dist = pdist2(xu , archive_xu);                                % this xu is upper level new infill xu, not added into archive_xu
[~, idx] = sort(dist);
% for i = 1: length(idx)
%     if size(lower_trg{idx(i)}, 1) > 500
%         close_id = i;
%         break
%     end
% end
close_id = idx(1);

close_optxu = archive_xu(close_id, :);
close_optxl = lower_xl(close_id, :);

if  global_halfflag   % use global half solutions, ditched
    global_xl         = lhsdesign(global_halfflagsize, prob.n_lvar, 'criterion','maximin','iterations',1000);
    global_xl         = repmat(prob.xl_bl, global_halfflagsize, 1) ...
        + repmat((prob.xl_bu - prob.xl_bl), global_halfflagsize, 1) .* global_xl;
    xu_ext            = repmat(xu, global_halfflagsize, 1);
    global_fl         = prob.evaluate_l(xu_ext, global_xl);
    size_localsamples = global_halfflagsize;
else
    global_xl = [];
    global_fl = [];
    size_localsamples = global_halfflagsize * 2;
end

close_xlarchvie = lower_trg{close_id}(:, 1:end-1);
[expensive_x, cheap_x, cheap_f, correlation, sigma] = sample_fromCloseLandscape(size_localsamples, prob, xu, close_optxu, close_optxl, close_xlarchvie,...
    'gaussian', true);

% [expensive_x, cheap_x, cheap_f, correlation] = sample_fromCloseLandscapeOnlyHalfOverlap(size_localsamples, prob, xu, close_optxu, close_optxl, close_xlarchvie,...
%    'gaussian', true);

end


function close_optxlcokring =  expandingHF(expensive_x, expensive_f, cheap_x, cheap_f, prob, xu, sigma, nc, close_optxu)
%
% 
% cheap solutions have already gone through repeat check
x_trg{1,:} = cheap_x;                                           % cheap
y_trg{1,:} = cheap_f;

% [1] identify center
c = 1;
[~, idx] = sort(expensive_f);
close_optxlcokring = expensive_x(idx(1), :);
cokrg_lb = min(expensive_x, [], 1);                              % use cheap because it expensive overlaps with cheap
cokrg_ub = max(expensive_x, [], 1);



while onboundary_check(close_optxlcokring, cokrg_lb, cokrg_ub, prob.xl_bl, prob.xl_bu)
    fprintf('HF expand x %d \n', c);
    
    [~, idx] = sort(expensive_f);
    current_bestx = expensive_x(idx(1), :);
    current_bestf = expensive_f(idx(1), :);
    
   
    xl_samples = [];
    while size(xl_samples, 1) < nc
        % [2] use neighbour parameters
        R = mvnrnd(current_bestx, sigma, nc * 2);                 % expand the same size as archive to avoid for loop
        
        mask_lbcomply = R > prob.xl_bl;
        mask_ubcomply = R < prob.xl_bu;
        mask_comply = mask_lbcomply .* mask_ubcomply; %
        mask_comply = sum(mask_comply, 2);
        
        idx = mask_comply == prob.n_lvar;
        
        tmp = R(idx, :);
        xl_samples = [xl_samples; tmp];

        
        % if boundary is not changing then no additional evaluation is
        % needed  
        tmp_xlexpand = [xl_samples; expensive_x];
        tmp_xllb = min(tmp_xlexpand, [], 1);
        tmp_xlub = max(tmp_xlexpand, [], 1);
        
        % whether boundary is changing ?
        % if the close_optxl is still on boundary after the expanding
        % samples, then this expanding is not needed
        if onboundary_checkTrh(close_optxlcokring, tmp_xllb, tmp_xlub, 1e-2, prob.xl_bl, prob.xl_bu)
            return
        end
  
        
    end
    
    % [3] repeat until close_optxlcokring stay within the bounds of HF
    xl_samples = xl_samples(1:nc, :);
    xuu = repmat(xu, nc, 1);
    fl_samples = prob.evaluate_l(xuu, xl_samples);
    expensive_x = [expensive_x; xl_samples];
    expensive_f = [expensive_f; fl_samples];
    
    
    [~, ia, ~]  = unique(expensive_x, 'rows');    % check unique
    expensive_x = expensive_x(ia, :);
    expensive_f = expensive_f(ia, :);
    fprintf('expensive candidate go through distance check \n');
    [expensive_x, expensive_f] = keepdistance_rawdata(expensive_x, expensive_f, prob.xl_bl, prob.xl_bu);
    
    x_trg{2,:} = expensive_x;
    y_trg{2,:} = expensive_f;
    
    
    % Fitting the cokring model
    co_mdl = oodacefit(x_trg, y_trg);
    
    cokrg_lb = min(x_trg{2,:} ,[], 1); % use cheap because it expensive overlaps with cheap
    cokrg_ub = max(x_trg{2,:} ,[], 1);
    [close_optxlcokring] = cokriging_decisionmaking(co_mdl, prob, cokrg_lb, cokrg_ub, 'hx',expensive_x, 'hf',expensive_f, 'lx',cheap_x ,'lf', cheap_f,...
        'boundary_updatesize',   floor(size(expensive_x, 1)/2),  'xu', xu, 'extend', false);
    
    % visualize
    visualise_expand(prob, xu, close_optxu, prob.xl_bl, prob.xl_bu,co_mdl,  expensive_x, expensive_f, xl_samples, close_optxlcokring, cheap_x, cheap_f);   
    c = c + 1;
    
    expensive_x = [expensive_x; close_optxlcokring]; 
    close_optflcokring = prob.evaluate_l(xu, close_optxlcokring);
    expensive_f = [expensive_f; close_optflcokring];
    
    
    [~, idx] = sort(expensive_f);
    close_optxlcokring = expensive_x(idx(1), :);
   
    
end


end
 



function visualise_expand(prob, xu, close_optxu, lb, ub, cokrg_mdl, co_expensivex, co_expensivey, new_samplex, cokrg_xl, co_cheapx, co_cheapy)
fignh               = figure(3);
fignh.Position      = [50, 100, 1000, 800];
nt                  = 100;


x1_tst              = linspace(lb(1), ub(1), nt);
x2_tst              = linspace(lb(2), ub(2), nt);
[msx1, msx2]        = meshgrid(x1_tst, x2_tst);
msx11 = msx1(:);
msx22 = msx2(:);
xl  = [msx11, msx22];
xum = repmat(xu, nt*nt, 1);
fl  = prob.evaluate_l(xum, xl);
fl  = reshape(fl, [nt, nt]);



% seed and returned xl in current landscape
subplot(2, 2, 1);
title('Current landscape');
surf(msx1, msx2, fl, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;
% seed_flnew = prob.evaluate_l(xu, close_optxl);
% scatter3(close_optxl(1), close_optxl(2), seed_flnew,  80, 'r', 'filled' ); hold on;

[seed_flcok, ~] = prob.evaluate_l(xu, cokrg_xl);
scatter3(cokrg_xl(:, 1), cokrg_xl(:, 2), seed_flcok,  40, 'blue', 'filled'); hold on;
scatter3(co_expensivex(:, 1), co_expensivex(:, 2), co_expensivey,  40, 'yellow', 'filled'); hold on;

legend('landscape', 'decision from cokrg', 'local search result from cokrg', 'expensive training', 'Location','north');



% seed in the orignal landscape
subplot(2, 2, 2);
title('close landscape');
xum = repmat(close_optxu, nt * nt, 1);
fl = prob.evaluate_l(xum, xl);
fl = reshape(fl, [nt, nt]);
% fl_seed = prob.evaluate_l(close_optxu, close_optxl);
surf(msx1, msx2, fl, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
% scatter3(close_optxl(1), close_optxl(2), fl_seed,  80, 'r', 'filled' ); hold on;
scatter3(co_cheapx(:, 1), co_cheapx(:, 2), co_cheapy,  40, 'yellow', 'filled'); hold on;
legend('landscape','cheap optimal', 'cheap train',   'Location','north');


% seed in surrogate landscape
subplot(2, 2, 3);
title('current landscape built from cokrg')
f = cokrg_mdl.predict(xl);
f = reshape(f, [nt, nt]);
surf(msx1, msx2, f, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;


[seed_fpred, ~] = cokrg_mdl.predict(cokrg_xl);
scatter3(cokrg_xl(1), cokrg_xl(2), seed_fpred,  40, 'b', 'filled' ); hold on;


% plot predicted results from expensive
expensive_predy = cokrg_mdl.predict(co_expensivex);
scatter3(co_expensivex(:, 1), co_expensivex(:, 2), expensive_predy,  80, 'r', 'filled' ); hold on;
[new_samplef, ~] = cokrg_mdl.predict(new_samplex);

scatter3(new_samplex(:, 1), new_samplex(:, 2), new_samplef, 40, 'g', 'filled'); hold on;

legend('cokrg landscape current', 'cokrg starting point', 'cok train', 'new samples', 'Location', 'north');

% seed in surrogate landscape
subplot(2, 2, 4);

% scatter3(close_archive(:, 1), close_archive(:, 2), close_archive(:, 3),  80, 'r', 'filled' ); hold on;
% pause(1);
close(fignh);
end

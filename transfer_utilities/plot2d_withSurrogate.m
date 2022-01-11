
function plot2d_withSurrogate(xu, lb, ub, match_xl, seed_xl, prob, seed_xu, delta_mdl, close_mdl, corr_xl)
% this function plot lower landscape 
% with starting point from surrogate landscape

fignh               = figure(3);
nt                  = 100;

cp                  = zeros(nt, nt);
x1_tst              = linspace(lb(1), ub(1), nt);
x2_tst              = linspace(lb(2), ub(2), nt);
[msx1, msx2]        = meshgrid(x1_tst, x2_tst);
msx11 = msx1(:);
msx22 = msx2(:);
xl = [msx11, msx22];
xum = repmat(xu, nt*nt, 1);
fl = prob.evaluate_l(xum, xl);
fl = reshape(fl, [nt, nt]);

fl_seed = prob.evaluate_l(xu, seed_xl);
fl_match = prob.evaluate_l(xu, match_xl);

% seed and returned xl in current landscape
subplot(2,2,1);
surf(msx1, msx2, fl, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;
scatter3(seed_xl(1), seed_xl(2), fl_seed,  80, 'r', 'filled' ); hold on;
scatter3(match_xl( 1), match_xl( 2), fl_match,  80, 'g', 'filled' ); hold on;
legend('landscape','seed', 'return', 'Location','north');

% seed in the orignal landscape
subplot(2, 2, 2);
xum = repmat(seed_xu, nt*nt, 1);
fl = prob.evaluate_l(xum, xl);
fl = reshape(fl, [nt, nt]);
fl_seed = prob.evaluate_l(seed_xu, seed_xl);
surf(msx1, msx2, fl, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;
scatter3(seed_xl(1), seed_xl(2), fl_seed,  80, 'r', 'filled' ); hold on;

fl_corr = prob.evaluate_l(seed_xu, corr_xl);
scatter3(corr_xl(1), corr_xl(2), fl_corr,  80, 'k', 'filled' ); hold on;

legend('landscape','seed in original', 'correlation seed in original', 'Location','north');

% seed in surrogate landscape
subplot(2, 2, 3);

param.GPR_type   = 2;
param.no_trials  = 1;

[delta_f, ~] = Predict_Objective(delta_mdl, xl, param);
[close_f, ~] = Predict_Objective(close_mdl, xl, param);

f = close_f + delta_f;
f =  reshape(f, [nt, nt]);
surf(msx1, msx2, f, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;

[delta_fseed, ~] = Predict_Objective(delta_mdl, seed_xl, param);
[close_fseed, ~] = Predict_Objective(close_mdl, seed_xl, param);
seed_fpred = close_fseed + delta_fseed;
scatter3(seed_xl(1), seed_xl(2), seed_fpred,  80, 'r', 'filled' ); hold on;
legend('surrogate landscape current','seed predicted', 'Location','north');

pause(1);
close(fignh);

end



% Predicting using a Kriging model
function [y_prd, s_prd] = Predict_Objective(gprMdl,x_tst,param)
y_prd=[]; s_prd=[];
for i=1:length(gprMdl)
    if(param.GPR_type==1)
        [tmp_mu,tmp_sigma] = predict(gprMdl{i},x_tst);
    else
        [tmp_mu,tmp_sigma] = Predict_DACE(gprMdl{i},x_tst,i);
    end
  
    y_prd = [y_prd tmp_mu];
    s_prd = [s_prd tmp_sigma];
end
end


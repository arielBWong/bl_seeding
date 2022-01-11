function plot2d_withNeighbourOptimal(xu, lb, ub, match_xl, close_optxl, prob, close_optxu)
% this function plot lower landscape 
% with starting point from surrogate landscape

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


fl_match = prob.evaluate_l(xu, match_xl);

% seed and returned xl in current landscape
subplot(2, 2, 1);
title('Current landscape');
surf(msx1, msx2, fl, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;
seed_flnew = prob.evaluate_l(xu, close_optxl);
scatter3(close_optxl(1), close_optxl(2), seed_flnew,  80, 'r', 'filled' ); hold on;
scatter3(match_xl(1), match_xl(2), fl_match,  80, 'g', 'filled' ); hold on;

legend('landscape', 'close optimal', 'found optimal', 'Location','north');



% seed in the orignal landscape
subplot(2, 2, 2);
title('close landscape');
xum = repmat(close_optxu, nt*nt, 1);
fl = prob.evaluate_l(xum, xl);
fl = reshape(fl, [nt, nt]);
fl_seed = prob.evaluate_l(close_optxu, close_optxl);
surf(msx1, msx2, fl, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;
seed_flneighbour = prob.evaluate_l(close_optxu, close_optxl);

scatter3(close_optxl(1), close_optxl(2), seed_flneighbour,  80, 'r', 'filled' ); hold on;

legend('neighbour landscape','seed in neighbour',   'Location','north');


pause(1);
close(fignh);

end


function [f] = current_landscapeF(x, mdl)

param.GPR_type   = 2;
param.no_trials  = 1;

[f, ~] = Predict_Objective(mdl, x, param);

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



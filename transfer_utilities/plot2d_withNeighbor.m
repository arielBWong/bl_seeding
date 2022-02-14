function plot2d_withNeighbor(xu, lb, ub, match_xl, close_optxl, prob, close_optxu)
% this function plot lower landscape 
% with starting point from surrogate landscape

fignh1 = figure(1);
% fignh.Position   = [50, 100, 1000, 400];
nt   = 100;


x1_tst  = linspace(lb(1), ub(1), nt);
x2_tst  = linspace(lb(2), ub(2), nt);
[msx1, msx2] = meshgrid(x1_tst, x2_tst);
msx11 = msx1(:);
msx22 = msx2(:);
xl  = [msx11, msx22];
xum = repmat(xu, nt*nt, 1);
fl  = prob.evaluate_l(xum, xl);
fl  = reshape(fl, [nt, nt]);


fl_match = prob.evaluate_l(xu, match_xl);

% seed and returned xl in current landscape

title('Current landscape');
surf(msx1, msx2, fl, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;
% 
flmatch = prob.evaluate_l(xu, match_xl);
s1 = scatter3(match_xl(1), match_xl(2), flmatch,  80, 'r', 'filled' ); hold on;

%
flstart = prob.evaluate_l(xu, close_optxl);
s2 = scatter3(close_optxl(1), close_optxl(2), flstart, 80, 'g', 'filled'); hold on;

lgd = legend([s1, s2],  'LL optimal solution',  'Neighbor optimal solution' , 'Location','northeast');
lgd.FontSize = 14;
lgd.FontWeight = 'bold';

xlabel('x_{l1}', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('x_{l2}', 'FontSize', 20, 'FontWeight', 'bold', 'rotation', 0);

hcb = colorbar;




pause(1);
close(fignh1);

% seed in the orignal landscape
fignh2 = figure(2);
title('close landscape');
xum = repmat(close_optxu, nt*nt, 1);
fl = prob.evaluate_l(xum, xl);
fl = reshape(fl, [nt, nt]);
fl_seed = prob.evaluate_l(close_optxu, close_optxl);
surf(msx1, msx2, fl, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;
%
flopt = prob.evaluate_l(close_optxu, close_optxl);
s1 = scatter3(close_optxl(1), close_optxl(2), flopt,  80, 'g', 'filled'); hold on;

lgd = legend([s1],   'Neighbor optimal solution' , 'Location','northeast');
lgd.FontSize = 14;
lgd.FontWeight = 'bold';

xlabel('x_{l1}', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('x_{l2}', 'FontSize', 20, 'FontWeight', 'bold', 'rotation', 0);
colorbar;
pause(1);
close(fignh2);

end




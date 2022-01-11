
function plot2dlower(xu, lb, ub, match_xl, seed_xl, prob, seed_xu)

fignh = figure(3);
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

subplot(1,2,1);
surf(msx1, msx2, fl, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;
scatter3(seed_xl(1), seed_xl(2), fl_seed,  80, 'r', 'filled' ); hold on;
scatter3(match_xl( 1), match_xl( 2), fl_match,  80, 'g', 'filled' ); hold on;
legend('landscape','seed', 'return', 'Location','north');

subplot(1,2,2);
xum = repmat(seed_xu, nt*nt, 1);
fl = prob.evaluate_l(xum, xl);
fl = reshape(fl, [nt, nt]);
fl_seed = prob.evaluate_l(seed_xu, seed_xl);
surf(msx1, msx2, fl, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;
scatter3(seed_xl(1), seed_xl(2), fl_seed,  80, 'r', 'filled' ); hold on;
legend('landscape','seed in original', 'Location','north');

pause(1);
close(fignh);


end



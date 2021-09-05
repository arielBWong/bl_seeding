
function seedinglower_plot2D(h1, h2, h3, prob, newxu,newxl, seed_xu, seed_xl, archive, lower_archive, idx)
cla(h1);
cla(h2);
cla(h3)

n = 200;
xu1 = linspace(prob.xu_bl(1), prob.xu_bu(1), n);
xu2 = linspace(prob.xu_bl(2), prob.xu_bu(2), n);
[xu11, xu22] = meshgrid(xu1, xu2);
xu11_all = xu11(:);
xu22_all = xu22(:);
xu_all = [xu11_all, xu22_all];
xl_all = prob.get_xlprime(xu_all);
fu_plot = prob.evaluate_u(xu_all, xl_all);
fu_plot = reshape(fu_plot, n, n);

axes(h1);
surf(xu11, xu22, fu_plot, 'FaceAlpha',0.5, 'EdgeColor', 'none' ); hold on;

fu_1 = prob.evaluate_u(newxu, newxl);
plot3(newxu(1), newxu(2), fu_1,   'Marker', 'o', 'MarkerSize', 8, 'Color', 'red','MarkerFaceColor' ,'red' ); hold on;
plot3(archive.xu(:, 1), archive.xu(:, 2), archive.fu,  'Marker', 'o', 'MarkerSize', 8, 'Color', 'yellow','MarkerFaceColor' ,'yellow'); hold on;

[seedfu, ~] = prob.evaluate_u(seed_xu, seed_xl);
plot3(seed_xu(1), seed_xu(2), seedfu,   'Marker', 'o', 'MarkerSize', 8, 'Color', 'green','MarkerFaceColor' ,'green' );


xl1 = linspace(prob.xl_bl(1), prob.xl_bu(1), n);
xl2 = linspace(prob.xl_bl(2), prob.xl_bu(2), n);
[xl11, xl22] = meshgrid(xl1, xl2);
xl11_all = xl11(:);
xl22_all = xl22(:);
xl_all = [xl11_all, xl22_all];
xu_all = repmat(newxu, n^2, 1);
fl_plot = prob.evaluate_l(xu_all, xl_all);
fl_plot = reshape(fl_plot, n, n);

axes(h2);

fl_1 = prob.evaluate_l(newxu, seed_xl);
fl_2 = prob.evaluate_l(newxu, newxl);

surf(xl11, xl22, fl_plot, 'FaceAlpha',0.5, 'EdgeColor', 'none' ); hold on;
plot3(seed_xl(1), seed_xl(2), fl_1,   'Marker', 'o', 'MarkerSize', 8, 'Color', 'green','MarkerFaceColor' ,'green'); hold on;
plot3(newxl(1), newxl(2), fl_2,   'Marker', 'o', 'MarkerSize', 8, 'Color', 'red','MarkerFaceColor' ,'red'); hold on;
plot3(lower_archive.init_xl(:,1), lower_archive.init_xl(:, 2), lower_archive.init_fl(end, :)','Marker', 'x');
r = corr(lower_archive.init_fl(end, :)', lower_archive.init_fl(idx(1), :)');
lower_archive.r = [lower_archive.r , r];

title(r);
hold on;

axes(h3);
xu_all = repmat(seed_xu, n^2, 1);
fl_plot = prob.evaluate_l(xu_all, xl_all);
fl_plot = reshape(fl_plot, n, n);
surf(xl11, xl22, fl_plot, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;

fl_3 = prob.evaluate_l(seed_xu, seed_xl);
plot3(seed_xl(1), seed_xl(2), fl_3,   'Marker', 'o', 'MarkerSize', 8, 'Color', 'green','MarkerFaceColor' ,'green');
plot3(lower_archive.init_xl(:,1), lower_archive.init_xl(:, 2), lower_archive.init_fl(idx(1), :)','Marker', 'x');
pause(1);
end
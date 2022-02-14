problem_folder = strcat(pwd,'/problems/TP3');
addpath(problem_folder);

prob = smd7mp(1, 1, 1);

fignh  = figure(1);
% fignh.Position   = [50, 100, 1000, 800];
nt = 100;
x1_tst   = linspace(prob.xu_bl(1), prob.xu_bu(1), nt);
x2_tst  = linspace(prob.xu_bl(2), prob.xu_bu(2), nt);
[msx1, msx2] = meshgrid(x1_tst, x2_tst);
msx11 = msx1(:);
msx22 = msx2(:);
xu  = [msx11, msx22];
xl = prob.get_xlprime(xu);
fu  = prob.evaluate_u(xu, xl);
fu = reshape(fu, [nt, nt]);

surf(msx1, msx2, fu, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;
s1 = scatter3(prob.xu_prime(1), prob.xu_prime(2), prob.fu_prime,  180, [0.9290 0.6940 0.1250], 'filled'); hold on;
xu = [0, 0.1];
xl_opt = prob.get_xlprime(xu);
fu2 = prob.evaluate_u(xu, xl_opt);
s2 = scatter3(xu(1), xu(2), fu2,  200, [0.8500 0.3250 0.0980],  'filled', 'd'); hold on;
lgd = legend([s1, s2],  'x_{u}^*', 'A Neighbor x_{u}', 'Location','north');
lgd.FontSize = 20;
lgd.FontWeight = 'bold';

xlabel('x_{u1}', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('x_{u2}', 'FontSize', 20, 'FontWeight', 'bold');
zlabel('F_{u}^*',  'FontSize', 20,  'FontWeight', 'bold', 'Rotation',0);

view(-70, 40);
pause(1);
close(fignh);




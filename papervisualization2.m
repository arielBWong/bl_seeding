problem_folder = strcat(pwd,'/problems/TP3');
addpath(problem_folder);

prob = smd7mp(1, 1, 1);
%----------------------------------------------------
fignh  = figure(1);
% fignh.Position   = [50, 100, 1000, 800];
nt = 100;

x1_tst   = linspace(prob.xl_bl(1), prob.xl_bu(1), nt);
x2_tst  = linspace(prob.xl_bl(2), prob.xl_bu(2), nt);
[msx1, msx2] = meshgrid(x1_tst, x2_tst);
msx11 = msx1(:);
msx22 = msx2(:);
xu = prob.xu_prime;
xuu = repmat(xu, nt*nt, 1);
xl = [msx11, msx22];

fl = prob.evaluate_l(xuu, xl);
fl = reshape(fl, [nt, nt]);
contour(msx1, msx2, fl); hold on;
s1 = scatter(prob.xl_prime(1), prob.xl_prime(2),  180, [0.9290 0.6940 0.1250], 'filled'); hold on;


xu2 = [0, 0.1];
xl2 = prob.get_xlprime(xu2);
fl = prob.evaluate_l(xu, xl2);
s2 =  scatter(xl2(1), xl2(2), 200, [0.8500 0.3250 0.0980], 'filled');

lgd = legend([s1, s2],  'x_{l}^* of x_u^*',  'x_{l}^{*} from x_{u}^{ne}' , 'Location','southeast');
lgd.FontSize = 20;
lgd.FontWeight = 'bold';
%legend boxoff

xlabel('x_{l1}', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('x_{l2}', 'FontSize', 20, 'FontWeight', 'bold', 'rotation', 0);




close(fignh);
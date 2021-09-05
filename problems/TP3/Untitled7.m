

n = 50;

xu = [1.5, 1.2, 0.5, 0.2, 0.1, 0.01, 0.001,];
prob= mp_bl();
xl1 = linspace(prob.xl_bl, prob.xl_bu, 200);

for i = 1:7
    cla;
    f = prob.evaluate_l(xu(i), xl1');
    plot(xl1, f,  'LineWidth', 2); hold on;
    xl = prob.get_xlprime(xu(i));
    f1 = prob.evaluate_l(xu(i), xl(1));
%     f2 = prob.evaluate_l(xu(i), xl(2));
    
    plot(xl(1), f1,    'Marker', 'o', 'MarkerSize', 8, 'Color', 'red','MarkerFaceColor' ,'red' );
    % plot(xl(2), f2,    'Marker', 'o', 'MarkerSize', 8, 'Color', 'red','MarkerFaceColor' ,'red' );
    
end


classdef smd6x
    properties
        p = 1;
        q = 1;
        r = 1;
        n_lvar;
        n_uvar;
        n_con;
        n_obj;
        xu_bl;
        xu_bu;
        xl_bl;
        xl_bu;
        ref;
        uopt = [];
        lopt = [];
        xu_prime = [];
        xl_prime = [];
        name = 'smd6x_levy';
    end
    methods
        function obj = smd6x(p, q, r)
            if nargin > 1
                obj.q = q;
                obj.p = p;
                obj.r = r;
            end
            obj.n_lvar = obj.q + obj.r;
            obj.n_uvar = obj.p + obj.r;

            % bounds
            %init bound upper level
            obj.xu_bl = zeros(1, obj.p + obj.r);
            obj.xu_bu = ones(1, obj.p + obj.r);
            
            
            % init bound lower level
            obj.xl_bl = [ones(1, obj.q) * -10, ones(1, obj.r)*-5];
            obj.xl_bu = ones(1, obj.q + obj.r) * 10 ;
            
            % best value
            obj.uopt     = obj.q ;
            obj.lopt     = 0;
            obj.xu_prime = [zeros(1, obj.p), zeros(1, obj.r)];
            obj.xl_prime = [ones(1, obj.q), zeros(1, obj.r)];
            

        end
        function [f, c] = evaluate_u(obj, xu, xl)
            xu1 = xu(:, 1 : obj.p);
            xu2 = xu(:, obj.p + 1: obj.p + obj.r);
            
            xl1 = xl(:, 1 : obj.q);
            xl2 = xl(:, obj.q + 1 : obj.q + obj.r);
            
            %-obj
            f = sum((xu1).^2, 2) ...   %F1
                + sum(xl1.^2, 2) ...   %F2 -- change to cooperative 
                + sum((xu2).^2, 2) ... %F3
                - sum((xu2 - xl2).^2, 2); % can be cooperative or conflict
            % -cie
            c = [];
            
        end
        
        function [f, c] = evaluate_l(obj, xu, xl)
           
            %------convergence difficulty levy
            x = xl(:, 1: obj.q, :);
            d = size(x, 2);
            
            for ii = 1:d
                w(:, ii) = 1 + (x(:, ii) - 1) ./ 4;
            end
            
            term1 = (sin(pi  .* w(:, 1))) .^2;
            term3 = (w(:, d)-1) .^ 2  .*  (1 + (sin(2 .* pi .*  w(:, d))) .^2);
            
            sum0 = 0;
            for ii = 1:(d-1)
                wi = w(:, ii);
                new = (wi - 1) .^ 2 .*  (1+10 .* (sin(pi .* wi + 1)) .^ 2);
                sum0 = sum0 + new;
            end
            
            f = term1 + sum0 + term3;
            
            %------------------------------------
            
            xu1 = xu(:, 1 : obj.p);
            xu2 = xu(:, obj.p + 1: obj.p + obj.r);
            
            xl1 = xl(:, 1 : obj.q);
            xl2 = xl(:, obj.q+1 : obj.q+obj.r);
            f1  = sum((xu1).^2, 2);                % function dependence 
            f3  = sum((xu2 - xl2).^2, 2);          % interaction
            f   = f1+ f + f3;
            c   = [];
            
        end
        
        
  
    
        
        function plot_problem(obj)
            m = 101;
            x1 = linspace(obj.xl_bl(1), obj.xl_bu(1), m);
            x2 = linspace(obj.xl_bl(2), obj.xl_bu(2), m);
            [x1, x2] = meshgrid(x1, x2);
            f = zeros(m, m);
            
            for i = 1:m
                for j = 1:m
                    f(i, j) = evaluate_l(obj, [], [x1(i, j), x2(i, j)]);
                end
            end
            surf(x1, x2, f); hold on;
            colormap jet
            shading interp
            
            
            
        end
    end
end
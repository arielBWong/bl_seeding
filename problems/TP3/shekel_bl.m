classdef shekel_bl
    properties
        n_lvar = 1;
        n_uvar = 1;
        m = 4;
        xu_bl;
        xu_bu;
        xl_bl;
        xl_bu;
        p;
        xprime;
        name = 'rastrigin';
        
    end
    methods
        function obj = shekel_bl()
            % level variables
            obj.n_lvar = 1;
            obj.n_uvar = 1;
            obj.p = pi/2 - 0.01; 
            
            % bounds
            %init bound upper level
            obj.xu_bl  = ones(1,obj.n_uvar) * obj.p * (-1);
            obj.xu_bu = ones(1, obj.n_uvar) * obj.p;
            
            % init bound lower level
            obj.xl_bl = ones(1, obj.n_lvar)  * -5;
            obj.xl_bu = ones(1, obj.n_lvar) * 6;
            
        end
        
        function [f, c] = evaluate_u(obj, xu, xl)
            %-obj
            f = sum(xu.^2, 2) + sum((xl-0.1*abs(tan(xu))).^2, 2);
            c = [];
        end
        
        
        function [f, c] = evaluate_l(obj, xu, xl)
            %-obj
              b = [0.05, 0.2, 0.3, 0.3]';
           C = [-4, 4, -2, 2;
                -4, 4, -2, 2;
                -4, 4, -2, 2;
                -4, 4, -2, 2];
            
            outer = 0;
            gap = 0;
            nv = 1;
            for ii = 1:obj.m
                bi = b(ii);
                inner = 0;
                for jj = 1:nv
                    xj = xl(:, jj);
                    g = (5-gap);
                    Cji = C(jj, ii) + g;
                    if Cji > 5 || Cji < -5
                        Cji = 0;
                    end
                    inner = inner + (xj-Cji).^2;
                end
                outer = outer + 1./(inner+bi);
            end
            outer = outer - sum(xl.^2, 2);
            f = -outer;
            %-con
            c = [];
        end
        
        
        function xl = get_xlprime(obj, xu)
            xl =  0.1*abs(tan(xu));
        end
    end
    
end

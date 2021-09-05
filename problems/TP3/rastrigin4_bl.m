classdef rastrigin4_bl
    properties
        n_lvar = 1;
        n_uvar = 1;
        xu_bl;
        xu_bu;
        xl_bl;
        xl_bu;
        p;
        xprime;
        name = 'rastrigin';
        
    end
    methods
        function obj = rastrigin4_bl()
            % level variables
            obj.n_lvar = 1;
            obj.n_uvar = 1;
            obj.p = pi/2 - 0.01; 
            
            % bounds
            %init bound upper level
            obj.xu_bl  = ones(1,obj.n_uvar) * obj.p * (-1);
            obj.xu_bu = ones(1, obj.n_uvar) * obj.p;
            
            % init bound lower level
            obj.xl_bl = ones(1, obj.n_lvar)  * 0;
            obj.xl_bu = ones(1, obj.n_lvar) * 10;
            
        end
        
        function [f, c] = evaluate_u(obj, xu, xl)
            %-obj
            f = sum((xu + obj.p).^2, 2) + sum((xl-0.1*abs(tan(xu))).^2, 2);
            c = [];
        end
        
        
        function [f, c] = evaluate_l(obj, xu, xl)
            %-obj
            f = 0;
            n = 1;
            for i = 1: obj.n_lvar
                k = xl(:, i) - 0.1 * abs(tan(xu(:, i)));
                fb = k.^2 - 20 * cos( pi * k.^2);
                f = f + fb;
            end
            
            f = 10 * obj.n_lvar + f;
            %-con
            c = [];
        end
        
        
        function xl = get_xlprime(obj, xu)
            xl =  0.1*abs(tan(xu));
        end
    end
    
end

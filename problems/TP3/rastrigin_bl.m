classdef rastrigin_bl
    properties
        n_lvar = 2;
        n_uvar = 2;
        xu_bl;
        xu_bu;
        xl_bl;
        xl_bu;
        xprime;
        name = 'rastrigin';

    end
    methods
        function obj = rastrigin_bl()
            % level variables
            obj.n_lvar = 2;
            obj.n_uvar = 2;
            
            % bounds
            %init bound upper level
            obj.xu_bl = zeros(1,obj.n_uvar);
            obj.xu_bu = ones(1, obj.n_uvar)*2;   
           
            % init bound lower level
            obj.xl_bl = ones(1, obj.n_lvar)  * -8;
            obj.xl_bu = ones(1, obj.n_lvar) * 8;   
            
        end
        
        function [f, c] = evaluate_u(obj, xu, xl) 
            %-obj
            f = sum(xu.^2, 2) + sum(xl.^2, 2);
            c = [];          
        end
        
        
        function [f, c] = evaluate_l(obj, xu, xl)         
            %-obj
            f = 0;
            n = 1;
           for i = 1: obj.n_lvar
               fb = (xl(:, i)- exp(xu(:, i))).^2 - 20 * cos(pi * (xl(:, i)-exp(xu(:, i))));
               f = f + fb;
           end
           
         f = 10 * obj.n_lvar + f;
            %-con
            c = [];
        end
       
        
        function xl = get_xlprime(obj, xu)
            xl = exp(xu);
        end
    end
    
end

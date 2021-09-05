classdef testproblem1
    properties
        n_lvar = 1;
        n_uvar = 1;
        xu_bl;
        xu_bu;
        xl_bl;
        xl_bu;
        uopt = 0;
        lopt = 0;
        name;
    end
    methods
        function obj = testproblem1()
          
            obj.name = 'testproblem1';
            
            % level variables
            obj.n_lvar = 1;
            obj.n_uvar = 1;
            
            % bounds
            % init bound upper level

            obj.xu_bl = -1;
            obj.xu_bu = 1;
           
            obj.xl_bl = -1;
            obj.xl_bu = 1;
            
        end
        
        function [f, c] = evaluate_u(obj, xu, xl)
            f = xu.^2 + xl.^2;
            c = [];
        end
        function [f, c] = evaluate_l(obj, xu, xl)
            f = (xu - xl).^2;
            %-cie
            c=[];
            
        end
        
        function xl_prime = get_xlprime(obj, xu)
          xl_prime = xu;
        end
    end
end

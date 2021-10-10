classdef smd7mp
    properties
        p = 1;
        q = 2;
        r  = 1;
        n_lvar;
        n_uvar;
        xu_bl;
        xu_bu;
        xl_bl;
        xl_bu;
        name;
        fu_prime;
        fl_prime ;
        xu_prime;
        xl_prime;

    end
    methods
        function obj = smd7mp(p, q, r)
            if nargin == 3
                obj.p = p;
                obj.q = q;
                obj.r = r;
            end
            obj.name = 'SMD7_mp';
            
            % level variables
            obj.n_lvar = obj.q + obj.r;
            obj.n_uvar = obj.p + obj.r;
            
            % bounds
            %init bound upper level
            xu_bl_1 = ones(1, obj.p)* (-5.0);
            xu_bu_1 = ones(1, obj.p) * 10.0;
            xu_bl_2 = ones(1, obj.r)* (-5.0);
            xu_bu_2 = ones(1, obj.r) * 1.0;
            obj.xu_bl = [xu_bl_1, xu_bl_2];
            obj.xu_bu = [xu_bu_1, xu_bu_2];
            
            % init bound lower level
            xl_bl_1 = ones(1, obj.q) * (-5.0);
            xl_bu_1 = ones(1, obj.q) * 10.0;
           %  xl_bl_2 = ones(1, obj.r) * (0 + 1e-6);
            % xl_bu_2 = ones(1, obj.r) * (exp(1));
            
            xl_bl_2 = ones(1, obj.r) * (0 );
            xl_bu_2 = ones(1, obj.r) * (10);
            
            obj.xl_bl = [xl_bl_1, xl_bl_2];
            obj.xl_bu = [xl_bu_1, xl_bu_2];
            
            obj.xu_prime = zeros(1, obj.n_uvar);
            obj.xl_prime = obj. get_xlprime(obj.xu_prime);
            obj.fu_prime = obj.evaluate_u(obj.xu_prime, obj.xl_prime);
            obj.fl_prime = obj.evaluate_l(obj.xu_prime, obj.xl_prime);
            
        end
        
        function [f, c] = evaluate_u(obj, xu, xl)
            xu1 = xu(:, 1 : obj.p);
            xu2 = xu(:, obj.p + 1: obj.p + obj.r);
            
            xl1 = xl(:, 1 : obj.q);
            xl2 = xl(:, obj.q+1 : obj.q+obj.r);
            
          
            xu_prime2 = zeros(1, obj.r);            
            reg = 5; % absolute value of  lower bound of xl2        
            t = 0.1 * tan(pi/2 - 0.01 - (abs(xu(:, obj.p+1 : end) )- xu_prime2) * pi/(2 * reg));
            
            %-obj
            m = [1:obj.p];
            
%             f = 1+1/400*sum((xu1).^2, 2) - prod(cos(xu1./sqrt(m)), 2) ...
%                 - sum((xl1).^2, 2) ...
%                 + sum((xu2).^2, 2) - sum((xu2 - log(xl2)).^2, 2);    

              f = 1+1/400*sum((xu1).^2, 2) - prod(cos(xu1./sqrt(m)), 2) ...
                  - sum((xl1).^2, 2) ...
                  + sum((xu2).^2, 2)...
                  - sum( (xl2 - t).^2, 2); 

              
            %-cie
            c = [];
        end
        
        function [f, c] = evaluate_l(obj, xu, xl)
            xu1 = xu(:, 1 : obj.p);
            xu2 = xu(:, obj.p + 1: obj.p + obj.r);
            
            xl1 = xl(:, 1 : obj.q);
            xl2 = xl(:, obj.q+1 : obj.q+obj.r);
          
            
            xu_prime2 = zeros(1, obj.r);            
            reg = 5; % absolute value of  lower bound of xl2        
            t = 0.1 * tan(pi/2 - 0.01 - (abs(xu(:, obj.p+1 : end) )- xu_prime2) * pi/(2 * reg));
            
            
            %-obj
%             f = sum((xu1).^3, 2) ...
%                     + sum((xl1).^2, 2) ...
%                     + sum((xu2 - log(xl2)).^2, 2);
            
            %dr = obj.xu_bu(obj.p + 1 :obj.p + obj.r) - obj.xu_bl(obj.p+1 :obj.p+obj.r);      
            f   = sum((xu1).^3, 2) ...
                    + sum((xl1).^2, 2) ...
                    + sum( (xl2 - t).^2, 2); 
                   
            % grafting moving peak problem
            fmp = mp_module(obj, xu, xl);
              
             f = f -  fmp;
            %-cie
            c = [];        
        end
        
        
        function xl_prime = get_xlprime(obj, xu)
            n = size(xu, 1);
            xl_prime = zeros(n, obj.n_lvar);
            for i = 1:obj.q
                xl_prime(:, i) = 0;
            end
            
           xu_prime2 = zeros(1, obj.r);            
           reg = 5; % absolute value of  lower bound of xl2        
           t = 0.1 * tan(pi/2 - 0.01 - (abs(xu(:, obj.p+1 : end) )- xu_prime2) * pi/(2 * reg));
            
            j = 1;
            for i = obj.q + 1 : obj.q + obj.r     
                % xl_prime(:,i) = exp(xu(:, obj.p+ j));             
                xl_prime(:, i) =t(:, j);
                      
                j = j + 1;
            end
        end
        
        function xl_local = get_otherlocal(obj, xu)
            xl_local = mp_othercenters(obj, xu);
        end
    end
end

classdef mp_bl
    properties
        n_lvar = 1;
        n_uvar = 1;
        xu_bl;
        xu_bu;
        xl_bl;
        xl_bu;
        p;
        xprime;
        name = 'mp_bl';
        
    end
    methods
        function obj = mp_bl()
            % level variables
            obj.n_lvar = 1;
            obj.n_uvar = 1;
            obj.p = pi/2 - 0.01; 
           
            
            % bounds
            %init bound upper level
            obj.xu_bl  = ones(1,obj.n_uvar) * obj.p * (0);
            obj.xu_bu = ones(1, obj.n_uvar) * obj.p;
            
            % init bound lower level
            obj.xl_bl = ones(1, obj.n_lvar)  * 0;
            obj.xl_bu = ones(1, obj.n_lvar) * 100;       
        end
        
        
        function [f, c] = evaluate_u(obj, xu, xl)
            %-obj
            f = sum(xu.^2, 2);
            c = [];
        end
        
        
        function [f, c] = evaluate_l(obj, xu, xl)
            %-obj
%             px = [8, 64, 67, 55, 4;
%                     20, 13, 76, 15, 7;
%                     9, 19, 27, 67, 24;
%                     66, 87, 65, 19, 43;
%                     76, 32, 43, 54, 65];
                              
                
               px = [66, 64, 67, 55, 4;
                        76, 13, 76, 15, 7;
                        55, 19, 27, 67, 24;
                        20, 87, 65, 19, 43;
                        8, 32, 43, 54, 65];
                    
              w = ones(5, 1) * 0.1;
              w(2) = 0.03;
              w(3) =0.05;
              h = ones(5, 1) * 50;
              h(2) = 30;
             h(3) = 20; 
              f=[];
              v = 0.6 * (tan(obj.p -xu));
              nx = size(xl, 1);
                
             v = 0.6* (tan(obj.p - xu));
             
              for i = 1:nx
                  fi = [];
                  for j = 1:3
                    fj =  h(j)/(1+  w(j) * sum((xl(i, :) - (px(j, 1:obj.n_lvar)  - v) ).^2, 2));
                    fi = [fi, fj];
                  end
                  f = [f; max(fi)];
              end             
            %-con
            c = [];
        end
        
        
        function xl = get_xlprime(obj, xu)
            
              px = [66, 64, 67, 55, 4;
                        76, 13, 76, 15, 7;
                        55, 19, 27, 67, 24;
                        20, 87, 65, 19, 43;
                        8, 32, 43, 54, 65];
             v = 0.9 *(tan(obj.p -xu));
             xl1 =  px(1, 1:obj.n_lvar)  - v;
             xl2 =  px(2, 1:obj.n_lvar)  - v;
            xl = xl1;
             
            % xl = [xl1, xl2];
        end
    end
    
end

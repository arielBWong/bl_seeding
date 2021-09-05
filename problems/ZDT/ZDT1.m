classdef ZDT1
    
    properties
        n_uvar;
        xu_bl;
        xu_bu;
        n_con;
        n_obj;
        name;
    end
    
    methods
        function obj = ZDT1(n_var)
            if nargin > 0
                obj.n_uvar = n_var;
                obj.xu_bl = zeros(1, n_var) ;
                obj.xu_bu = ones(1, n_var);
                obj.n_con = 0;
                obj.n_obj = 2;
                obj.name = 'ZDT1';
            else
                obj.n_uvar = 6;
                obj.xu_bl = zeros(1, obj.n_var) ;
                obj.xu_bu = ones(1, obj.n_var);
                obj.n_con = 0;
                obj.n_obj = 2;
                obj.name = 'ZDT1';
            end
        end
        
        function [y, c] = evaluate_u(obj, PopDec)
            % function only for test, incomplete implementation
            PopObj(:,1) = PopDec(:,1);
            g = 1 + 9*mean(PopDec(:,2:end),2);
            h = 1 - (PopObj(:,1)./g).^0.5;
            PopObj(:,2) = g.*h;
            y = PopObj;
            c = [];
        end
        
        function y = upper_pf(obj, N)
            R(:,1) = linspace(0,1,N)';
            R(:,2) = 1 - R(:,1).^0.5;
            y = R;
        end
    end
end
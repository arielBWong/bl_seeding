classdef bleaq_tp3
    properties
        p = 1;
        q = 0;
        r = 1;
        n_lvar;
        n_uvar;
        xu_bl;
        xu_bu;
        xl_bl;
        xl_bu;
        uopt = 0;
        lopt = 0;
        name;
        fu_prime = -18.6787;
        fl_prime = -1.0156;
    end
    methods
        function obj = bleaq_tp3()
            obj.name = 'bleaq_tp3';
            
            % level variables
            obj.n_lvar = 2;
            obj.n_uvar = 2;
            
            % bounds
            %init bound upper level
            obj.xu_bl =  zeros(1,obj.n_uvar);  
            obj.xu_bu = 10*ones(1,obj.n_uvar);   
            
            % init bound lower level   
            obj.xl_bl = zeros(1, obj.n_lvar); 
            obj.xl_bu = 10*ones(1,obj.n_lvar);   
            
        end
        
        function [f, c] = evaluate_u(obj, xu, xl)
            %At optima x1=0,x2=2,y1=1.8750,xl2=0.9062,fu=18.6787,fl=1.0156
             x1 = xu(:, 1);
             x2 = xu(:, 2);
             y1 = xl(:, 1);
             y2 = xl(:, 2);
             functionValue = -x1.^2-3*x2.^2-4*y1+y2.^2;
        
             f = functionValue;
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Write the constraints here
             equalityConstrVals = [];
             inequalityConstrVals(:, 1) = x1.^2+2*x2-4;
            %Lower level constraints included at upper level
            inequalityConstrVals(:, 2) = -3-x1.^2+2*x1-x2.^2+2*y1-y2;
            inequalityConstrVals(:, 3) = 4-x2-3*y1+4*y2;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            c=inequalityConstrVals;
            
        end
        function [f, c] = evaluate_l(obj, xu, xl)
            x1 = xu(:, 1);
            x2 = xu(:, 2);
            y1 = xl(:, 1);
            y2 = xl(:, 2);
            functionValue = 2*x1.^2+y1.^2-5*y2;

             f = functionValue;

             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Write the constraints here
            equalityConstrVals = [];
            inequalityConstrVals(:, 1) = -3-x1.^2+2*x1-x2.^2+2*y1-y2;
            inequalityConstrVals(:, 2) = 4-x2-3*y1+4*y2;
            c=inequalityConstrVals;            
        end
        
      
    end
end

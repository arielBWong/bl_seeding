function [xsol,fval,history, output] = lowerlevel_fmincon(x0, maxFE, lb, ub,funh_obj, funh_con)
% Set up shared variables with outfun
% Should also save the referenced 
history.x = [];
history.fval = [];

options = optimoptions(@fmincon,'OutputFcn',@outfun,...
    'Display','off',...
    'MaxFunctionEvaluations', maxFE);

[xsol,fval, ~, output] = fmincon(funh_obj,x0,[],[],[],[],lb,ub,funh_con,options);

    function stop = outfun(x,optimValues,state)
        
        stop = false;
        switch state
            case 'iter'
                % Concatenate current point and objective function
                % value with history. x must be a row vector.
                history.fval = [history.fval; optimValues.fval];
                history.x = [history.x; x];
            otherwise
        end
        
    end

end


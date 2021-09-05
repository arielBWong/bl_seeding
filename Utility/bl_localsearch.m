% upper level needs a local search 
function [xu, xl, n_up, n_low] = bl_localsearch(prob, xu_start, nest_parameter)
%---------------------------------


% to avoid repeated calculation of [f, c]
global xu_g
global xl_g
global fu_g
global cu_g
global ll_n
xu_g =[];
xl_g =[];
fu_g =[];
cu_g =[];
ll_n = 0;


% select a starting point and do local search
fmin_obj = @(x)blobj(x, prob, nest_parameter);
fmin_con = @(x)blcon(x, prob, nest_parameter);

opts                = optimset('fmincon');
opts.Algorithm      = 'sqp';
opts.Display        = 'off';
opts.MaxFunctionEvaluations ...
                    = nest_parameter.maxUpperFE;
[xu, ~, ~, output]  = fmincon(fmin_obj, xu_start, [], [],[], [],  ...
    prob.xu_bl, prob.xu_bu, fmin_con,opts);

[xl,~, ~]           = check_exist(xu, prob, nest_parameter.arc_con);
if isempty(xl)
   error('bl_local search has no xl stored');
end    
n_up                = output.funcCount;
n_low               = ll_n;
end


function fu =  blobj(xu, prob, nest_parameter)
global xu_g
global xl_g
global fu_g
global cu_g
global ll_n

[xl, fu, ~] = check_exist(xu, prob, nest_parameter.arc_con);
if ~isempty(xl)
   return 
end    
[xl, n]             = llmatch_trueEval(xu, prob, nest_parameter);
[fu, cu]            = prob.evaluate_u(xu, xl);

xu_g                = [xu_g; xu];
xl_g                = [xl_g; xl];
fu_g                = [fu_g; fu];
cu_g                = [cu_g; cu];
ll_n                = ll_n + n;

end

function [cu, cequ] =  blcon(xu, prob, nest_parameter)
cequ = [];

global xu_g
global xl_g
global fu_g
global cu_g
global ll_n

[xl, ~, cu] = check_exist(xu, prob, nest_parameter.arc_con);
if ~isempty(xl)
   return 
end
[xl, n]             = llmatch_trueEval(xu, prob, nest_parameter);

[fu, cu]            = prob.evaluate_u(xu, xl);

xu_g                = [xu_g; xu];
xl_g                = [xl_g; xl];
fu_g                = [fu_g; fu];
cu_g                = [cu_g; cu];
ll_n                = ll_n + n;

end

function [xl, fu, cu] = check_exist(xu, prob, constraint)
xl = []; fu = []; cu = [];


global xu_g
global xl_g
global fu_g
global cu_g

if isempty(xu_g)
    return;
else    
    diff = xu_g - xu; 
    ind = sum(diff, 2) == 0;
    
    if sum(ind)>0           % should be only 1, if exists
        % fprintf('found xu in global save %d\n', sum(ind));
        
        if sum(ind)>1
            name = prob.name;
            fprintf('there should not be repeated x more than once, prob, %s \n', name);
            nn = 1:size(xu_g, 1);
            archive_xug = nn(ind);
            ind = archive_xug(1);
        end
        xl = xl_g(ind, :);
        fu = fu_g(ind, :);
        if ~isempty(constraint)
            cu = cu_g(ind, :); else
            cu = [];
        end
    end
end
end

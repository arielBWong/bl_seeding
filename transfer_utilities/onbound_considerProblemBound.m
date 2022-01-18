function [flag] = onbound_considerProblemBound(x, lb, ub, prob)
% this function  check in normalized space
% whether x is on boundary

flag_up = false;
flag_down = true;

x_norm = (x - prob.xl_bl)./(prob.xl_bu - prob.xl_bl);
lb_norm = (lb - prob.xl_bl) ./(prob.xl_bu - prob.xl_bl);
ub_norm = (ub - prob.xl_bl) ./(prob.xl_bu - prob.xl_bl);

prob_unorm = ones(1, prob.n_lvar);
prob_lnorm = zeros(1, prob.n_lvar);

check_ub = ub_norm - x_norm;
id_ub = check_ub < 1e-6;
if any(id_ub > 0) 
   check_probound = prob_unorm - x_norm;
   if any(check_probound < 1e-6)
       flag_up = false; % if the second time it is still on problem bound, then consider it not on bound
   else
       flag_up = true; 
   end
 
end

check_lb = x_norm - lb_norm;
id_lb = check_lb < 1e-6;
if any(id_lb > 0)
    check_probound = x_norm - prob_lnorm;
    if any(check_probound < 1e-6)
        flag_up = false; % if the second time it is still on problem bound, then consider it not on bound
    else
        flag_up = true;
    end
end

flag = flag_up|flag_down;

end 
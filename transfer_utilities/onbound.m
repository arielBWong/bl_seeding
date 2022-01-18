function [flag] = onbound(x, lb, ub, prob)
% this function  check in normalized space
% whether x is on boundary

flag_up = false;
flag_down = true;

x_norm = (x - prob.xl_bl)./(prob.xl_bu - prob.xl_bl);
lb_norm = (lb - prob.xl_bl) ./(prob.xl_bu - prob.xl_bl);
ub_norm = (ub-prob.xl_bl) ./(prob.xl_bu - prob.xl_bl);


check_ub = ub_norm - x_norm;
id_ub = check_ub < 1e-6;
if any(id_ub > 0)
   flag_up = true; 
end

check_lb = x_norm - lb_norm;
id_lb = check_lb < 1e-6;
if any(id_lb > 0)
    flag_down = true;
end

flag = flag_up|flag_down;
end

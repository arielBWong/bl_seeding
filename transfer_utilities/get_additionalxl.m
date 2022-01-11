function [top_xl, top_fl] = get_additionalxl(prob, close_xlarchvie, seed_xu, close_xlprime)

[n, m] = size(close_xlarchvie);
seed_xus = repmat(seed_xu, n, 1);

close_fl = prob.evaluate_l(seed_xus, close_xlarchvie);

% start to locate a local region
prim     =  bump_detection(close_xlarchvie, close_fl);
nboxes   =  length(prim.boxes);

% plot the figure of each box
if m == 2               % only for 2d
   fig = figure(1); 
   scatter(close_xlarchvie(:, 1), close_xlarchvie(:, 2),  40, 'yellow', 'filled'); hold on;  % 2d plot
   
   lb = prob.xl_bl;
   ub = prob.xl_bu;
   nt = 100;
   x1_tst  = linspace(lb(1), ub(1), nt);
   x2_tst = linspace(lb(2), ub(2), nt);
   [msx1, msx2] = meshgrid(x1_tst, x2_tst);
   msx11 = msx1(:);
   msx22 = msx2(:);
   xl  = [msx11, msx22];
   xum = repmat(seed_xu, nt*nt, 1);
   fl = prob.evaluate_l(xum, xl);
   fl = reshape(fl, [nt, nt]);
   contour(msx1, msx2, fl); hold on;
   scatter(close_xlprime(1), close_xlprime(2), 80, 'red', 'filled' ); hold on;
end

% 
upper_stack = [];
lower_stack = [];
for ii = 1: nboxes
    [bump_lb, bump_ub] = boxboundary(prim, ii, prob);
    upper_stack        = [upper_stack; bump_ub];
    lower_stack        = [lower_stack; bump_lb];
    
    if m == 2
        % plot the current box
        rectangle('position', [bump_lb(1), bump_lb(2), bump_ub(1), bump_ub(2)], 'LineWidth',3); hold on;
    end
    
end

if m == 2 
    pause(1);
    close(fig);
end


top_xl = [];
top_fl = [];
end



function[bump_lb, bump_ub] = boxboundary(prim, ii, prob)
mentioned_var = prim.boxes{ii}.vars;
bump_lb = prob.xl_bl;
bump_ub = prob.xl_bu;

nb = length(mentioned_var);
if ~isempty(mentioned_var)
    for jj = 1:nb % varible indicator
        if ~isnan(prim.boxes{ii}.min(jj))
            bump_lb(mentioned_var(jj)) =  prim.boxes{ii}.min(jj);
        end
        
        if ~isnan(prim.boxes{ii}.max(jj))
            bump_ub(mentioned_var(jj)) =  prim.boxes{ii}.max(jj);
        end 
    end
end
end

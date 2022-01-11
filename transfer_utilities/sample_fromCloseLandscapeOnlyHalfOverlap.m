function[expensive_x, cheap_x, cheap_f, r] = sample_fromCloseLandscapeOnlyHalfOverlap(ns, prob, xu, seed_xu, seed_xl, close_xlarchive, varargin)
% This function samples from cheap landscape providing local information for
% future cokriging 
% 
p = inputParser;
addRequired(p, 'ns');                                        % number of expensive samples
addRequired(p, 'prob');
addRequired(p, 'xu');
addRequired(p, 'seed_xu');
addRequired(p, 'seed_xl');
addRequired(p, 'close_xlarchive');
addParameter(p, 'box', false);
addParameter(p, 'box_currentx', []);
addParameter(p, 'gaussian', true);                           % probing archive
addParameter(p, 'visual', false);
parse(p, ns, prob, xu, seed_xu, seed_xl, close_xlarchive, varargin{:});
% -----------------------------------
ns = p.Results.ns;  
prob = p.Results.prob;
xu = p.Results.xu;
seed_xu = p.Results.seed_xu;
seed_xl = p.Results.seed_xl;
close_xlarchive = p.Results.close_xlarchive;
box = p.Results.box;
box_currentx = p.Results.box_currentx;
gaussian = p.Results.gaussian;
visual = p.Results.visual;
% ==============

cm = 40; % cheap ones
[n, m] = size(close_xlarchive);
seed_xus = repmat(seed_xu, n, 1);
close_fl = prob.evaluate_l(seed_xus, close_xlarchive);  % just for shortcut, closefl can be stored

if n < ns
    fprintf('Problem %s, close archive number smaller than required number %d, shrink HF size\n', prob.name, ns);
    ns = n;
end

if box && gaussian
    error('two option cannot be true at the same time')
end

if ~box && ~gaussian
    error('two option cannot be false at the same time')
end

if box % no more in use
    corners = [seed_xl; box_currentx];
    corner_lb = min(corners, [], 1);
    corner_ub = max(corners, [], 1);
    
    xl_samples        = lhsdesign(ns, prob.n_lvar, 'criterion','maximin','iterations',1000);
    xl_samples        = repmat(corner_lb, ns, 1) ...
        + repmat((corner_ub - corner_lb), ns, 1) .* xl_samples; 
end



if gaussian
    sigma = covariation(close_xlarchive, seed_xl);
    xl_samples = [];
    nc = size(close_xlarchive, 1);
    
    while size(xl_samples, 1) < ns 
        if nc < ns
            error('lf has less than HF');
        end
        R = mvnrnd(seed_xl, sigma, nc);                 % expand the same size as archive to avoid for loop
        
        mask_lbcomply = R > prob.xl_bl;
        mask_ubcomply = R < prob.xl_bu;
        mask_comply   = mask_lbcomply .* mask_ubcomply; % 
        mask_comply   = sum(mask_comply, 2); 
        
        idx = mask_comply == prob.n_lvar;
        
        xl_samples = R(idx, :);
        
    end

end

% k=1 smallest pairwise distance to observations in close_xlarchive for each observation in xl_samples
% [~, I] = pdist2(close_xlarchive, xl_samples,  'euclidean', 'Smallest', 1);

expensive_x = [];
cheap_x = [];
cheap_f = [];
for i = 1:floor(ns/2)
    [~, I] = pdist2(close_xlarchive, xl_samples(i, :),  'euclidean', 'Smallest', 1); % process one by one
    
    expensive_x = [expensive_x; close_xlarchive(I, :)];
    cheap_x = [cheap_x; close_xlarchive(I, :)];
    cheap_f = [cheap_f; close_fl(I, :)];
      
    % unique process, to make sure use as much archive as possible
    close_xlarchive(I, :) = [];
    close_fl(I, :) = [];
end

for i = floor(ns/2)+1: ns
    [~, I] = pdist2(close_xlarchive, xl_samples(i, :),  'euclidean', 'Smallest', 1); % process one by one
    
    expensive_x = [expensive_x; xl_samples(i, :)];
    cheap_x = [cheap_x; close_xlarchive(I, :)];
    cheap_f = [cheap_f; close_fl(I, :)];
      
    % unique process, to make sure use as much archive as possible
    close_xlarchive(I, :) = [];
    close_fl(I, :) = [];
end


xuu = repmat(xu, ns, 1);
expensive_f = prob.evaluate_l(xuu, expensive_x); % lazy step

r = corr(cheap_f(1: floor(ns)), expensive_f(1:floor(ns)));
fprintf('Correlation with closest landscape is %0.4f \n', r);
% in case close_xlarchive is global variable
% copy it out

if nc - ns < cm
    fprintf('Problem %s, close archive size after HF extraction smaller than required number %d, shrink LF (eliminated those used for HF) size is ', prob.name, cm)
    cm = nc - ns;
    fprintf('%d \n', cm);
end


cheap_x = [cheap_x; close_xlarchive(1: cm, :)];
cheap_f = [cheap_f; close_fl(1:cm, :)];

if visual
    visualization(prob, m, close_xlarchive, seed_xl, seed_xu, expensive_x);
end
end

function flag = inbound(x, lb, ub)
flag = true;
if any(x - lb < 0) || any(ub - x < 0)
    flag = false;
end 
end


function[] = visualization(prob, m, close_xlarchive, close_xlprime, seed_xu, expensive_x)
if m == 2               % only for 2d
    fig = figure(1);
    scatter(close_xlarchive(:, 1), close_xlarchive(:, 2),  80, 'yellow', 'filled'); hold on;  % 2d plot
    
    lb = prob.xl_bl;
    ub = prob.xl_bu;
    nt = 100;
    x1_tst = linspace(lb(1), ub(1), nt);
    x2_tst = linspace(lb(2), ub(2), nt);
    [msx1, msx2] = meshgrid(x1_tst, x2_tst);
    msx11 = msx1(:);
    msx22 = msx2(:);
    xl  = [msx11, msx22];
    xum = repmat(seed_xu, nt*nt, 1);
    fl = prob.evaluate_l(xum, xl);
    fl = reshape(fl, [nt, nt]);
    contour(msx1, msx2, fl); hold on;
    scatter(close_xlprime(1), close_xlprime(2), 40, 'red', 'filled' ); hold on;
    scatter(expensive_x(:, 1), expensive_x(:, 2), 80, 'Marker', 'X', 'LineWidth', 2);
    
    
    
    legend('cheap samples', 'landscape', 'cleap optimal', 'selection');
    pause(1);
    close(fig);
    
else
    error('This visualization does not serve more than 2d');
end
end
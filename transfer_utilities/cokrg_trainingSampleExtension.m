function [expensive_x, cheap_x, cheap_f, correlation, close_optxu, close_optxl, global_xl, global_fl, sigma, close_id] = cokrg_trainingSampleExtension(xu , ...
    archive_xu, lower_trg, prob, samplesize, lower_xl)
% this function identify close landscape for kriging (the one with global search)
% and identify training data for both Hx and Lx
dist = pdist2(xu , archive_xu);                                % this xu is upper level new infill xu, not added into archive_xu
[~, idx] = sort(dist);
for i = 1: length(idx)
    if size(lower_trg{idx(i)}, 1) > 500
        close_id = idx(i);
        break
    end
end
% close_id = idx(1);

global_xl = [];
global_fl = [];

seed_xu = archive_xu(close_id, :);
seed_xl = lower_xl(close_id, :);

close_optxu = seed_xu;
close_optxl = seed_xl;

ns = samplesize * 2;
close_xlarchive = lower_trg{close_id}(:, 1:end-1);


sigma = covariation(close_xlarchive, seed_xl);
nc = size(close_xlarchive, 1);
std_2 = sum( (close_xlarchive - seed_xl) .* (close_xlarchive - seed_xl), 1 )/ (nc - 1);
std_radius = sqrt(sum(std_2));
xl_samplesIn = [];
xl_samplesOut = [];
xl_samples = [];

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
    tmp = R(idx, :);
    
    % only add samples outside standard deviation
    xl_samplesDelta = abs(tmp - seed_xl);
    radius_Delta = sqrt(sum(xl_samplesDelta .* xl_samplesDelta, 2));
    
    mask_stdOut = radius_Delta > std_radius; 
    mask_stdIn = radius_Delta < std_radius;
    
    xl_samplesIn = tmp(mask_stdIn, :);
    xl_samplesOut = tmp(mask_stdOut, :);
    
    if size(xl_samplesIn, 1) < samplesize || size(xl_samplesOut, 1) < samplesize
        error('problem in sampling');
    end
    
    xl_samples = [xl_samplesIn(1:samplesize, :); xl_samplesOut(1:samplesize, :)];
    
end

close_xuu = repmat(seed_xu, nc, 1);
close_fl = prob.evaluate_l(close_xuu, close_xlarchive);

expensive_x = [];
cheap_x = [];
cheap_f = [];
for i = 1:ns
    [~, I] = pdist2(close_xlarchive, xl_samples(i, :),  'euclidean', 'Smallest', 1); % process one by one
    
    expensive_x = [expensive_x; close_xlarchive(I, :)];
    cheap_x = [cheap_x; close_xlarchive(I, :)];
    cheap_f = [cheap_f; close_fl(I, :)];
    
    % unique process, to make sure use as much archive as possible
    close_xlarchive(I, :) = [];
    close_fl(I, :) = [];
end

xuu = repmat(xu, ns, 1);
expensive_f = prob.evaluate_l(xuu, expensive_x); % lazy step

correlation = corr(cheap_f, expensive_f);
fprintf('Correlation with closest landscape is %0.4f \n', correlation);

end



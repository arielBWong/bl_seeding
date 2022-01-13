function [expensive_x, expensive_f, cheap_x, cheap_f, correlation] = cokrg_trainingSampleExtension(xu, prob, ...
    archive_xu, archive_xl, lower_searchdata_cly, samplesize)
% this function identify close landscape for kriging (the one with global search)
% and identify training data for both Hx and Lx
% input:
% xu - current upper variable
% prob - bl problem
% archive_xu - upper searched xu  
% archive_xl - matching lower xl
% lower_searchdata_cly - each xu's full lower search data
% samplesize - lower level search size
% output:
% expensive_x: HF data samples
% cheap_x: LF data samples
% cheap_f: LF objective values 
% correlation: correlation between LF and HF function values
%--------------------------------------------------------------------------

[~,~, close_id] = retrieve_neighbour(xu, lower_searchdata_cly, archive_xu, archive_xl);
ns = samplesize;

% identifying LF points (all evaluated)
close_xlarchive = lower_searchdata_cly{close_id}(:, 1:end-1);
seed_xl = archive_xl(close_id, :);

sigma = covariation(close_xlarchive, seed_xl);
nc = size(close_xlarchive, 1);
std_2 = sum( (close_xlarchive - seed_xl) .* (close_xlarchive - seed_xl), 1 )/ (nc - 1);  % sum over each element
std_radius = sqrt(sum(std_2));
xl_samplesIn = [];
xl_samplesOut = [];
xl_samples = [];

while size(xl_samples, 1) < ns
    if nc < ns
        error('Lf has less than HF');
    end
    R = mvnrnd(seed_xl, sigma, nc);                 % expand the same size as archive to avoid for loop
    
    mask_lbcomply = R > prob.xl_bl;
    mask_ubcomply = R < prob.xl_bu;
    mask_comply = mask_lbcomply .* mask_ubcomply; %
    mask_comply = sum(mask_comply, 2);
    
    idx = mask_comply == prob.n_lvar;
    tmp = R(idx, :);
    
    % only add samples outside standard deviation
    xl_samplesDelta = abs(tmp - seed_xl);
    radius_Delta = sqrt(sum(xl_samplesDelta .* xl_samplesDelta, 2));
    
    mask_stdOut = radius_Delta > std_radius; 
    mask_stdIn = radius_Delta < std_radius;
    
    xl_samplesIn = tmp(mask_stdIn, :);
    xl_samplesOut = tmp(mask_stdOut, :);
    
    if size(xl_samplesIn, 1) < int16(samplesize/2 )|| size(xl_samplesOut, 1) <int16( samplesize/2)
        % because take half within std, half outside std 
        error('problem in sampling, valid sample size is smaller than need');
    end
    
    xl_samples = [xl_samplesIn(1:int16(samplesize/2), :); xl_samplesOut(1:int16(samplesize/2), :)];
    
end

% close_xuu = repmat(seed_xu, nc, 1);  
% close_fl = prob.evaluate_l(close_xuu, close_xlarchive);
close_fl = lower_searchdata_cly{close_id}(:,  end);

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
expensive_f = prob.evaluate_l(xuu, expensive_x); 

correlation = corr(cheap_f, expensive_f);
fprintf('Correlation with closest landscape is %0.4f \n', correlation);

end



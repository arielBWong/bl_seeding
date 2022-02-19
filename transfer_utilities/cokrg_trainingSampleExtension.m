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
seed_xu = archive_xu(close_id, :);

sigma = covariation(close_xlarchive, seed_xl);
nc = size(close_xlarchive, 1);
nc = 1000;

% https://au.mathworks.com/help/matlab/ref/std.html
std_2 = sum( (close_xlarchive - seed_xl) .* (close_xlarchive - seed_xl), 1 )/ (nc - 1);  % sum over each element
std_2 = sqrt(std_2);

xl_samplesIn = [];
xl_samplesOut = [];
xl_samples = [];

bound_count = 0;

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

    mask2d_stdIn =  (abs(tmp - seed_xl) < std_2); 
    mask_stdIn = sum(mask2d_stdIn, 2) == prob.n_lvar;
    mask_stdOut = ~mask_stdIn;
    
    xl_samplesIn = tmp(mask_stdIn, :);
    xl_samplesOut = tmp(mask_stdOut, :);
    
    if size(xl_samplesIn, 1) < int16(samplesize/2 )|| size(xl_samplesOut, 1) <int16( samplesize/2)
        % because take half within std, half outside std
        % error('problem in sampling, valid sample size is smaller than need');
        fprintf('sample size is too small \n');
        if size(xl_samplesIn, 1) < int16(samplesize/2 )
            xl_samples = [xl_samples;  xl_samplesIn];
        else
            xl_samples = [xl_samples; xl_samplesIn(1:int16(samplesize/2), :)];
        end
        %------------
        if size(xl_samplesOut, 1) <int16( samplesize/2)
            xl_samples = [xl_samples; xl_samplesOut ];
        else
             xl_samples = [xl_samples;  xl_samplesOut(1:int16(samplesize/2), :) ];
        end
        
        % recount sample size
        ns = size(xl_samples, 1);
        
        
    else
         xl_samples = [xl_samplesIn(1:int16(samplesize/2), :); xl_samplesOut(1:int16(samplesize/2), :)];
    end
    
    % xl_samples = [xl_samplesIn(1:int16(samplesize/2), :); xl_samplesOut(1:int16(samplesize/2), :)];
    % check on which dimension  with std will violate boundary
    

    check_upperbound = (seed_xl + std_2) >  prob.xl_bu;
    check_lowerbound = (seed_xl - std_2) < prob.xl_bl;

    if  any(check_upperbound) 
        new_boundsample =  seed_xl .* (~check_upperbound) + prob.xl_bu .* (check_upperbound);
        % xl_samples(int16(samplesize/2) + 1, :) = []; % delete a middle point
        
        rnd_idx = randperm(ns);
        xl_samples(rnd_idx(1), :) = []; % delete a random point
        
        xl_samples = [xl_samples; new_boundsample];
        bound_count = bound_count + 1;
    end

    if any(check_lowerbound) 
        new_boundsample = seed_xl .* (~check_lowerbound) + prob.xl_bl .* check_lowerbound;
        % xl_samples(int16(samplesize/2) + 1, :) = []; % delete a middle point
        
        rnd_idx = randperm(ns);
        xl_samples(rnd_idx(1), :) = [];   % delete a random point
         
        
        xl_samples = [xl_samples; new_boundsample];
        bound_count = bound_count + 1;
    end   
end

close_fl = lower_searchdata_cly{close_id}(:,  end);

expensive_x = [];
cheap_x = [];
cheap_f = [];
for i = 1:ns-bound_count
    [~, I] = pdist2(close_xlarchive, xl_samples(i, :),  'euclidean', 'Smallest', 1); % process one by one
    
    expensive_x = [expensive_x; close_xlarchive(I, :)];
    cheap_x = [cheap_x; close_xlarchive(I, :)];
    cheap_f = [cheap_f; close_fl(I, :)];
    
    % unique process, to make sure use as much archive as possible
    close_xlarchive(I, :) = [];
    close_fl(I, :) = [];
end

for i = 1:bound_count
    expensive_x = [expensive_x; xl_samples(end - i+1, :)];
end

% create visualization, 
% visualize_std(seed_xl, std_2, prob, xu, expensive_x,  lower_searchdata_cly{close_id}(:, 1:end-1) );    


xuu = repmat(xu, ns, 1);
expensive_f = prob.evaluate_l(xuu, expensive_x); 

correlation = corr(cheap_f, expensive_f(1:ns-bound_count,:));   % correlation does not consider boundary
fprintf('Correlation with closest landscape is %0.4f \n', correlation);

% add more LF samples
extra_lb = min(expensive_x, [], 1);                         % use cheap because it expensive overlaps with cheap
extra_ub = max(expensive_x, [], 1);
cheap_extrax = lhsdesign(ns, prob.n_lvar, 'criterion','maximin','iterations',1000);
cheap_extrax = repmat(extra_lb, ns, 1) ...
    + repmat((extra_ub -extra_lb), ns, 1) .* cheap_extrax;

% match extra with archive 
for i = 1:ns
    [~, I] = pdist2(close_xlarchive, cheap_extrax(i, :),  'euclidean', 'Smallest', 1); % process one by one

    cheap_x = [cheap_x; close_xlarchive(I, :)];
    cheap_f = [cheap_f; close_fl(I, :)];
    
    % unique process, to make sure use as much archive as possible
    close_xlarchive(I, :) = [];
    close_fl(I, :) = [];
end


% plotsamples(seed_xu, seed_xl, prob, prob.xl_bl, prob.xl_bu,  lower_searchdata_cly{close_id}(:, 1:end-1), xl_samples, expensive_x, cheap_x);


end

function plotsamples(seed_xu, seed_xl, prob, lb, ub, archive_xl, samples_xl, HFx, LFx)
fignh = figure(1);
nt = 100;

cp = zeros(nt, nt);
x1_tst = linspace(lb(1), ub(1), nt);
x2_tst = linspace(lb(2), ub(2), nt);
[msx1, msx2] = meshgrid(x1_tst, x2_tst);
msx11 = msx1(:);
msx22 = msx2(:);
xl = [msx11, msx22];
xum = repmat(seed_xu, nt*nt, 1);
fl = prob.evaluate_l(xum, xl);
fl = reshape(fl, [nt, nt]);

xum = repmat(seed_xu, size(archive_xl, 1), 1);
fl_archive = prob.evaluate_l(xum, archive_xl);

xum = repmat(seed_xu, size(samples_xl, 1), 1);
fl_samples = prob.evaluate_l(xum, samples_xl);

xum = repmat(seed_xu, size(HFx, 1), 1);
fl_HFx = prob.evaluate_l(xum, HFx);

xum = repmat(seed_xu, size(LFx, 1), 1);
fl_LFx = prob.evaluate_l(xum, LFx);



scatter3(HFx(:, 1), HFx(:, 2), fl_HFx, 40, 'g', 'filled' ); hold on;

scatter3(LFx(:, 1), LFx(:, 2), fl_LFx, 70, 'k' ); hold on;
surf(msx1, msx2, fl, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
scatter3(archive_xl(:, 1), archive_xl(:, 2), fl_archive, 20, 'r', 'filled' ); hold on;

pause(1);


end

function visualize_std(mean, std, prob, xu, samplex, archivex)

fighn = figure(1);

nt = 100;
tst1 = linspace(prob.xl_bl(1), prob.xl_bu(1), nt);
tst2 = linspace(prob.xl_bl(2), prob.xl_bu(2), nt);

[tst11, tst22] = meshgrid(tst1, tst2);

tst11L = tst11(:);
tst22L = tst22(:);

xuu = repmat(xu, size(tst11L, 1), 1);
truef = prob.evaluate_l(xuu, [tst11L, tst22L]);
truef = reshape(truef, [nt, nt]);

subplot(1,2,1);
[M,c] = contour(tst11, tst22, truef); hold on;
c.LineWidth = 3;
scatter(mean(1), mean(2),  180, [0.9290 0.6940 0.1250], 'filled'); hold on;

lower_point = mean - std;
rectangle('Position', [lower_point, std(1)*2, std(2)*2]);

scatter(archivex(:, 1), archivex(:, 2), 200, 'g', 'filled'); hold on;
scatter(samplex(:, 1), samplex(:, 2), 120, 'r', 'filled'); hold on;

subplot(1, 2, 2);
surf(tst11, tst22, truef, 'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on;

pause(0.5);
close(fighn);

end



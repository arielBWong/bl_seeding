
function [localx, localf, localc, krg, krgc, arc_obj, arc_c] = klocalpoints(x, f, c, varargin)
% create a good local region around feasible min f
% not compatible with MO problem
% return local regions' data, model built from training data
%----------------------------
[best_x, ~, ~, ~, ~] = localsolver_startselection(x,  f,  c); %compatible with c
numx                 = size(x, 2);

currentk             = 3 * (numx + 1);
maxk                 = size(x, 1);

% initial train-test collection
idx                  = knnsearch(x, best_x, 'K', currentk);
localx               = x(idx, :);
localf               = f(idx, :);
if ~isempty(c)
    localc           = c(idx, :); else
    localc           = [];
end

% pick up test data/ and test data stay the same for the rest processes
trainx_num      = floor(0.8 * currentk);
p               = randperm(currentk);
test_idx        = p(trainx_num + 1 :end);
tstx            = localx(test_idx, :);
tstf            = localf(test_idx, :);

if ~isempty(c)
    tstc        = localc(test_idx, :); else
    tstc        = [];
end

% check trainx's theta compatible
% check trainx's uniqueness
[localx, localf, localc, train_idx, currentk] ...
                           = data_prepare(currentk, trainx_num, x, f, c, best_x,...
                                                localx, localf, localc, tstx, tstf, tstc);

% surrgoate model build, normalization using all training and test
[krg, krgc,arc_obj, arc_c] = surrogate_build(localx, localf, localc, train_idx);    
% predict on test data
current_mse                = surrogate_valid(krg, krgc, tstx, tstf, tstc, localf, arc_obj, arc_c);
previous_mse               = inf;


breakflag                  = true;
while breakflag
        
    if isnan(current_mse) % exception handling
        krg  = {};
        krgc = {};
        return;
    end
 
    if current_mse < previous_mse         
        if currentk == maxk
            break
        end
        
        % save previous output   
        previous_mse               = current_mse;
        previous_localx            = localx;
        previous_localf            = localf;
        previous_localc            = localc;
        previous_krg               = krg;
        previous_krgc              = krgc;
        previous_arcOBJ            = arc_obj;
        previous_arcCON            = arc_c;
        
        
        % expand k nearest neighbour
        currentk                    = currentk + 1;
        idx                         = knnsearch(x, best_x, 'K', currentk);
        localx                      = x(idx, :);
        localf                      = f(idx, :);
        if ~isempty(c)
            localc                  = c(idx, :); else
            localc                  = [];
        end
        % examine and repair training data quality
        [localx, localf, localc, train_idx, currentk] ...
                                    = data_prepare(currentk, trainx_num, x, f, c, best_x,...
                                                localx, localf, localc, tstx, tstf, tstc);
                                            
        % surrgoate model rebuild, (normalization using localf (train + test))
        [krg, krgc,arc_obj, arc_c]  = surrogate_build(localx, localf, localc, train_idx);    
        % predict on test data
        current_mse                 = surrogate_valid(krg, krgc, tstx, tstf, tstc, localf, arc_obj, arc_c);
        
    else
        breakflag                   = false;
        localx                      = previous_localx;
        localf                      = previous_localf;
        localc                      = previous_localc;
        krg                         = previous_krg;
        krgc                        = previous_krgc;
        arc_obj                     = previous_arcOBJ;
        arc_c                       = previous_arcCON;  
    end  
end

end


function [localx, localf, localc, train_idx, currentk] = data_prepare(currentk, trainx_num, x, f, c, best_x,...
                                                localx, localf, localc, tstx, tstf, tstc)
% from given local data, extract training data 
% make sure trainx's theta compatible to dace
% make sure trainx's uniqueness
[maxk, numx]            = size(x);
 
if size(localx, 1) ~= currentk
    error('k nearest neighbour, k is not same as localx/f/c size');
end

idx_tst                 = my_ismember(tstx, localx);
train_idx               = 1:currentk;
train_idx(idx_tst)      = [];
trainx                  = localx(train_idx, :);
trainf                  = localf(train_idx, :);
if ~isempty(c)
    trainc              = localc(train_idx, :); else
    trainc              = [];
end


notcomply_theta         = true;
while notcomply_theta && currentk < maxk
    
    % first check uniqueness
    if ~repeatpoint_check(trainx)
        [trainx, ia, ~] = unique(trainx, 'rows');
        trainf          = trainf(ia, :);
        if ~isempty(c)
            trainc      = trainc(ia, :); else
            trainc      = [];
        end
        % if fails unique check
        % local x f c excludes repeat points
        localx          = [trainx; tstx];
        localf          = [trainf; tstf];
        localc          = [trainc; tstc];
        train_idx       = 1: size(trainx, 1);
        
    end
    
  
    if size(trainx, 1) < trainx_num  % train size cannot be too small after unique process
        currentk         = currentk + 1;
        
        if currentk > maxk
            error('should not happen, training data can always cover ub < 1e-6 condition ');
        end
        idx              = knnsearch(x, best_x, 'K', currentk);
        localx           = x(idx, :);
        localf           = f(idx, :);
        if ~isempty(c)
            localc       = c(idx, :); else
            localc       = [];
        end
        % exclude test data
        idtst            = my_ismember(tstx, localx);
        train_idx        = 1:currentk;
        train_idx(idtst) = [];
        % pick up train data and test uniqueness
        trainx           = localx(train_idx, :);
        trainf           = localf(train_idx, :);
        if ~isempty(c)
            trainc       = localc(train_idx, :); else
            trainc       = [];
        end
        
    else
        % if currentk > (2 * (numx + 1))
        %     fprintf('local data comply with theta optimization and unique check with %d local samples\n', currentk);
        % end
        notcomply_theta  = false;
    end
end
end


function unique_check = repeatpoint_check(X)
% check training data uniqueness
% return processed dataset
% ------------------
[m, n] = size(X);
% Calculate distances D between points
mzmax = m*(m-1) / 2;        % number of non-zero distances
ij = zeros(mzmax, 2);       % initialize matrix with indices
D = zeros(mzmax, n);        % initialize matrix with distances
ll = 0;
for k = 1 : m-1
  ll = ll(end) + (1 : m-k);
  ij(ll,:) = [repmat(k, m-k, 1) (k+1 : m)']; % indices for sparse matrix
  D(ll,:) = repmat(X(k,:), m-k, 1) - X(k+1:m,:); % differences between points
end
if  min(sum(abs(D),2) ) == 0
  fprintf('Multiple design sites are detected in training data, training data needs adjustment \n'); 
  unique_check = false;
else
  unique_check = true;
end    
end



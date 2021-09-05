
function[localx, localf, localc, train_idx] = data_prepare(localx, localf, localc)
% check data quality
% uniqueness and theta upper bound
flag = repeatpoint_check(localx);
if ~flag
   
    [localx, ia, ~] = unique(localx, 'rows');
    localf          = localf(ia, :);
    if ~isempty(localc)
        localc      = localc(ia, :); 
    end
end

nx                  = size(localx, 1); % in case size changes
train_idx           = 1:nx;

end

function [unique_check] = repeatpoint_check(X)
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
  D(ll,:)   = repmat(X(k,:), m-k, 1) - X(k+1:m,:); % differences between points
end
if  min(sum(abs(D),2) ) == 0
  unique_check = false;
else
  unique_check = true;
end    
end
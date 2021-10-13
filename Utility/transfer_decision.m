function switch_indicator = transfer_decision(xu, archive, lower_archive, train_fl)
switch_indicator = struct();
dist = pdist2(xu, archive.xu);      % this xu is upper level new infill xu, not added into archive
[~, idx] = sort(dist);
r = corr(train_fl, lower_archive.init_fl(idx(1), :)');
switch_indicator.pearson = r;

% distance of Y
switch_indicator.dist = sqrt(sum((train_fl - lower_archive.init_fl(idx(1), :)').^2));

% skewness of Y
switch_indicator.skewness = [skewness(train_fl), skewness(lower_archive.init_fl(idx(1), :)')] ;

% kurtosis of Y
switch_indicator.kurtosis = [kurtosis(train_fl), kurtosis(lower_archive.init_fl(idx(1), :)')];

% length scale maybe


end


function [trgx, trgf] = keepdistance_rawdata(trgx, trgf, lbx, ubx)


n = size(trgx, 1);

normx = (trgx - lbx)./ (ubx - lbx);
keeplist = 1:n;   % this list is used to inidicate which ids are kept in terms of distance

d2 = pdist2(normx, normx);

% set lower triangle as big as possible
d2(logical(tril(ones(n, n)))) = inf;
L = d2 < 1e-3;
L = any(L, 1); % if one column has nonzero means the vector represented by this column ID should be eliminated

keeplist(L) = [];



% use keeplist as selection criterion
trgx = trgx(keeplist, :);
trgf = trgf(keeplist, :);

if length(keeplist) < n
    fprintf('too close distance discovered in preparing cokring \n ');
    fprintf('training data left size %d \n', size(trgx, 1));
end
end
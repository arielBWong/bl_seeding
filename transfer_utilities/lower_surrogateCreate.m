function [mdl] = lower_surrogateCreate(trgx, trgf, lbx, ubx)


n = size(trgx, 1);

normx = (trgx - lbx)./ (ubx - lbx);
keeplist = 1:n;   % this list is used to inidicate which ids are kept in terms of distance

d2 = pdist2(normx, normx);

% set lower triangle as big as possible
d2(logical(tril(ones(n, n)))) = inf;
L = d2 < 1e-6;
L = any(L, 1); % if one column has nonzero means the vector represented by this column ID should be eliminated

keeplist(L) = [];



% use keeplist as selection criterion
trgx = trgx(keeplist, :);
trgf = trgf(keeplist, :);


% train model 
% create  mapping from xu to fl
param.GPR_type      = 2;
param.no_trials     = 1;
mdl                 = Train_GPR(trgx, trgf, param);

end






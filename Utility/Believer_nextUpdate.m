function[best_x, info] = Believer_nextUpdate(train_x, train_y, xu_bound, xl_bound, ...
    num_pop, num_gen, train_c, normhn, varargin)
%----------------------------------
if ~isempty(varargin)
    fighn    = varargin{1};
    prob     = varargin{1};
end

if size(train_y, 2) > 1
    train_y_norm = normhn(train_y);
else
    train_y_norm = normalization_z(train_y);
end




param.GPR_type   = 2;
param.no_trials  = 1;
krg_obj          = Train_GPR(train_x, train_y_norm, param);

if ~isempty(train_c)
    krg_con      = Train_GPR(train_x, train_c, param); else
    krg_con      = {};
end

arc_obj.x        = train_x;
arc_obj.muf      = train_y_norm;

arc_con.x        = train_x;
arc_con.muf      = train_c;

% KB objective function
funh_obj         = @(x)Predict_GPR(krg_obj, x, param, arc_obj);
if ~isempty(train_c)
    funh_con     = @(x)Predict_GPR(krg_con, x, param, arc_con); else
    funh_con     = @(x)NoConstraint(x);
end

param_ea.gen     = num_gen;
param_ea.popsize = num_pop;

l_nvar           = size(train_x, 2);
[~,~,~, archive] = gsolver(funh_obj, l_nvar,  xl_bound, xu_bound, [], funh_con, param_ea);

newx             = archive.pop_last.X(1, :);

newx             = boundary_check(newx, xu_bound, xl_bound);
[best_x]         = surrogate_localsearch(newx, xu_bound, xl_bound, krg_obj, krg_con, arc_obj, arc_con);

best_x           = boundary_check(best_x, xu_bound, xl_bound);


%-----
info                = struct();
info.krg            = krg_obj;
info.arc_obj        = arc_obj;     
% info.stable = true;
if  ~isempty(train_c)
    info.krgc       = krg_con;
    info.arc_c      = arc_con ;
else
    info.krgc       = {};
    info.arc_c      = [];
end
%----
end


function f = NoConstraint(x)
f = [];
end

function[newx] = surrogate_localsearch(newx, xu_bound, xl_bound, krg_obj, krg_con, arc_obj, arc_con)
%----------------------------

if length(krg_obj) > 1
    % no local search for MO
    fprintf('no local surrogate search for MO \n');
     return
end

param.GPR_type   = 2;
param.no_trials  = 1;
% KB objective function
funh_obj          = @(x)Predict_GPR(krg_obj, x, param, arc_obj);
funh_con         = @(x)Predict_GPR(krg_con, x, param, arc_con);



opts                        = optimset('fmincon');
opts.Algorithm        = 'sqp';
opts.Display            = 'off';
opts.MaxFunctionEvaluations = 100;
[newx_local, ~, ~, ~] ...
                            = fmincon(funh_obj, newx, [], [],[], [],  ...
                                xl_bound,  xu_bound, funh_con, opts);
                            
% boundary check
newx_local                  = boundary_check(newx_local, xu_bound, xl_bound);
newx_localf                 = Predict_GPR(krg_obj, newx_local, param, arc_obj);
newx_localc                 = Predict_GPR(krg_con, newx_local, param, arc_con);                           

newf                        = Predict_GPR(krg_obj, newx, param, arc_obj);
newc                        = Predict_GPR(krg_con, newx, param, arc_con);

% decide which x to use
xx                         = [newx; newx_local];
ff                         = [newf; newx_localf];
cc                         = [newc; newx_localc];
[newx, ~, ~, ~, ~]         = localsolver_startselection(xx,  ff, cc);

end



function[localx, localf, localc, train_idx, continue_flag] = data_prepare(localx, localf, localc)
% check data quality
% uniqueness and theta upper bound
[flag, remove_index] = repeatpoint_check(localx);
if ~flag
    fprintf('before unique x size: %d \n', size(localx, 1));
    localx(remove_index, : ) = [];
    localf(remove_index, : ) = [];
    % [localx, ia, ~] = unique(localx, 'rows');
    % localf          = localf(ia, :);
    
    fprintf('after unique x size: %d \n', size(localx, 1));
    if ~isempty(localc)
        % localc      = localc(ia, :); 
        localc(remove_index, : ) = []; else
        localc      = [];
    end
end

nx                  = size(localx, 1);
train_idx           = 1:nx;

range               = minmax(localx');
ub                  = range(:,2);
if any(ub < 1e-6)
    continue_flag   = true; % theta range violate, outer process need expand local data size
    error('training data for determine dace theta [lb, ub] fails');
else
    continue_flag   = false;
    
end
end

function [unique_check, index] = repeatpoint_check(X)
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
  fprintf('Multiple design sites are detected in normalization \n'); 
  unique_check = false;
  repeat_ij= sum(abs(D),2) == 0;
  check_repeats = ij(repeat_ij, :); 
  
  index = [];
  for i = 1: size(check_repeats, 1)
      index = [index, check_repeats(i, 2)];
  end
  index = unique(index);
  fprintf('remove data number %d \n',  length(index));
else
  unique_check = true;
  index = [];
end    
end










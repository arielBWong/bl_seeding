
function [krg, krgc, arc_obj, arc_con] = surrogate_build(x, f, c, train_idx)
% this surrogate build considers training and test
% normalization with all test+train f
% but train only on train

if size(f, 2) > 1
    error('MO surrogate built not implemented');
end


% param.GPR_type=1 for GPR of Matlab; 2 for DACE
param.GPR_type   = 2;
param.no_trials  = 1;

f_norm      = normalization_z(f);

trainx      = x(train_idx, :);
trainy_norm = f_norm(train_idx, :);


arc_obj.x    = trainx;
arc_obj.muf  = trainy_norm;


if ~isempty(c)

    trainc      = c(train_idx, :);
    krgc        = Train_GPR(trainx, trainc, param);
    arc_con.x   = trainx;
    arc_con.muf = c(train_idx, :);
else
    krgc        = {};
    arc_con.x   = [];
    arc_con.muf = [];
end
krg             = Train_GPR(trainx, trainy_norm, param);
end



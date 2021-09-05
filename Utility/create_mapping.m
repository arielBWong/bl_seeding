

function [krg_obj, param] = create_mapping(xu, fl)
% upper lower mapping
if size(fl, 2) > 1
    fl_norm = normhn(fl);
else
    fl_norm = normalization_z(fl);
end

% create  mapping from xu to fl
param.GPR_type   = 2;
param.no_trials     = 1;
krg_obj                 = Train_GPR(xu, fl_norm, param);

end
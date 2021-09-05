%-----auxiliary function ---
function [f, s] = surrogate_predict(x, mdl, arc)


% param.GPR_type=1 for GPR of Matlab; 2 for DACE
param.GPR_type   = 2;
param.no_trials  = 3;

[f, s]= Predict_GPR(mdl, x, param, arc);
%  if mdl is empty, f=[], s=[]
end


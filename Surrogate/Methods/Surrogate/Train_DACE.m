function dmodel = Train_DACE(X,Y,no_trials)
%ref [1]: P16/ MATLAB Parametric Empirical Kriging (MPErK) User’s Guide

% Check design points
[m1, nx] = size(X);
[m2, ~]  = size(Y);
if m1 ~= m2
	error('X and Y must have the same number of rows')
end
range = minmax(X');
lb    = max(1e-6,range(:,1));
ub    = range(:,2);

Mi    = zeros(1, nx);
mi    = zeros(1, nx);


for i    = 1:nx  % ref [1]
    colx = X(:, i);
    Mi(i) = max(pdist(colx));
    tmp  = pdist(colx);
    tmp(tmp==0) ...
         = [];
    if size(tmp, 2) == 0
        % disp(colx);
        Mi(i) = - log(0.99^(1/nx)) / 1e-6;
        mi(i) = - log(0.01^(1/nx)) / 1e3  ;
    else
        mi(i) = min(tmp);
    end
    % fprintf('%f, ', (min(tmp)));
    
end

lb       = - log(0.99^(1/nx))./Mi; lb  = lb';
ub       = - log(0.01^(1/nx))./mi; ub  = ub';

fmin=Inf;
for trials = 1:no_trials
    theta_start  = lb+(ub-lb).*rand(nx,1);
    dmodel       = Fit_DACE(X, Y, @Regpoly0, @Corrgauss, theta_start', lb, ub);
    samples      = size(X, 1);
    %----- prediction
    if samples > 1
        [f, ssqr]    = Predictor(X, dmodel);
    else
        [f, ~, ssqr] = Predictor(X, dmodel);
    end
    if size(f,1)     == size(X,1)
        y = f;
    else
        y = f';
    end
    sigma = sqrt(abs(ssqr));
    y_prd = y;
    y_sig = sigma;
    %------------------
    
    fpr   = sqrt(sum((y_prd - Y).^2));    
    if(fpr<=fmin)
        output_model=dmodel;
        fmin=fpr;
    end
end

dmodel = output_model;
return
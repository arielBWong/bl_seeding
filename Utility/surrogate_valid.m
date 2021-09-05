

function current_mse = surrogate_valid(krg, krgc, tstx, tstf, tstc, localf, arc_obj, arc_con)
% check accuracy on test data
% single objective 
% training f normalization used train+test data, so denormalization also
% use all 
if size(localf, 2)>1
    error('validation not built for for MO');
end
[tstx, ia, ~]   = unique(tstx, 'row');
tstf            = tstf(ia, :);
[tst_pred, ~]   = surrogate_predict(tstx, krg, arc_obj);
tst_pred        = denormzscore(localf, tst_pred);
diff            = abs(tstf - tst_pred);


if ~isempty(krgc)
    [con_pred, ~]   = surrogate_predict(tstx, krgc, arc_con);
    diff_con        = abs(con_pred - tstc);
    
    num_c           = length(krgc);
    for i = 1: num_c        
        diff        = [diff; diff_con(:, i)];
    end
end

current_mse = sqrt(sum(diff.^2));
end

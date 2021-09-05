function train_y_norm = normalization_y(train_y_norm)
[num_y, d] = size(train_y_norm);

% check max min same mulfunctioning
if  d>1
    minf  = min(train_y_norm, [], 1);
    maxf = max(train_y_norm, [], 1);
    if num_y > 1 &&  all( (maxf - minf) >0) % 
        train_y_norm = (train_y_norm - minf) ./ (maxf - minf);
    else
         train_y_norm = train_y_norm;
    end 
    return 
end

% avoid same value on all ys
if num_y > 1 && (max(train_y_norm)-min(train_y_norm))>0
    train_y_norm = (train_y_norm -repmat(min(train_y_norm),num_y,1))./...
        repmat(max(train_y_norm)-min(train_y_norm),num_y,1);
else
    train_y_norm = train_y_norm;
end
end
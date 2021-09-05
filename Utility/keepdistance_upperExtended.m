
function [x, f, c, lx,  lf,  lc, flags] = keepdistance_upperExtended(x, f, c, lx,  lf,  lc, flags, ubx, lbx)
% eliminate x that is too close
% normalize first

normx               =  (x - lbx)./ (ubx - lbx);
newx                 =  normx(end, :);
oldx                   = normx(1:end-1, :); 

% check distance between x
d2                  = pdist2(newx, oldx);
[dmin, did]      = min(d2);

if dmin<1e-6          % process due to too close
    selectx             = [newx; normx(did, :)];
    selectf             = [f(end, :); f(did, :)];
    if isempty(c)       
        selectc         = [];
    else
        selectc         = [c(end, :); c(did, :)];
    end
    
    if size(selectx, 1) ~= 2
        error('compared x should be only 2 ');
    end
    [~, ~, ~ , ~, did2] ...
                        = localsolver_startselection(selectx, selectf, selectc);
                    
    if did2 == 1 % select new point
        % adjust archive 
        x(did, :) = []; f(did, :) = [];  lx(did, :) = [];   lf(did, :) = []; flags(did, :) = [];
        if ~isempty(c)
            c(did, :) = []; 
        end
        if ~isempty(lc)
            lc(did, :) = []; 
        end     
    else % select old point
         % adjust archive 
        x(end, :) = []; f(end, :) = [];  lx(end, :) = [];   lf(end, :) = []; flags(end, :) = [];
        if ~isempty(c)
            c(end, :) = []; 
        end
        if ~isempty(lc)
            lc(end, :) = []; 
        end     
    end         
end

end
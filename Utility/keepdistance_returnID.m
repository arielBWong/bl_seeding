
function [id] = keepdistance_returnID(x, f, c, ubx, lbx)
% id <0, no elimination
% id > 0; eliminate
% eliminate x that is too close
% normalize first



normx               =  (x - lbx)./ (ubx - lbx);
newx                =  normx(end, :);
newf                =  f(end, :);

oldx                =  normx(1:end-1, :);
oldx_orig           =  x(1:end-1, :);
oldf                =  f(1:end-1, :);

if ~isempty(c)
    oldc            = c(1:end-1, :);else
    oldc            = [];
end

d2                  = pdist2(newx, oldx);
[dmin, did]         = min(d2);

if dmin<1e-6 % process
    selectx             = [newx; normx(did, :)];
    selectf             = [newf; f(did, :)];
    if isempty(c)       
        selectc         = [];
    else
        selectc         = [c(end-1, :); c(did, :)];
        oldc(did, :)    = [];
    end
    
    oldx(did, :)        = [];
    oldx_orig(did, :)   = [];
    oldf(did, :)        = [];
    
    if size(selectx, 1) ~= 2
        error('compared x should be only 2 ');
    end
    [~, newf, newc, ~, did2] ...
                        = localsolver_startselection(selectx, selectf, selectc);
                    
    if did2 == 1
        newx = x(end, :);
        fprintf('distance fail, keep new point\n');
        id = did; % old to delete
    else
        newx = x(did, :);
        fprintf('distance fail, keep old point\n');
        id = size(x, 1); % new to delete
    end
    
    
          
else
      id = -1;
end

end

function [xk, xk_extra, fk, ck] = keepdistance_upper(x, f, c, x_extra, ubx, lbx)
% eliminate x that is too close
% normalize first


normx               =  (x - lbx)./ (ubx - lbx);
newx                 =  normx(end, :);
newf                 =  f(end, :);

oldx                   =  normx(1:end-1, :);
oldx_orig           =  x(1:end-1, :);
oldxe_orig          =  x_extra(1:end-1, :);
oldf                   =  f(1:end-1, :);

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
    oldxe_orig(did, :)  = [];
    oldf(did, :)        = [];
    
    if size(selectx, 1) ~= 2
        error('compared x should be only 2 ');
    end
    [~, newf, newc, ~, did2] ...
                        = localsolver_startselection(selectx, selectf, selectc);
                    
    if did2 == 1
        newx  = x(end, :);
        newxe = x_extra(end, :);
        fprintf('distance fail, keep new point\n');
    else
        newx  = x(did, :);
        newxe = x_extra(did, :);
        fprintf('distance fail, keep old point\n');
    end
    
    xk = [oldx_orig; newx]; fk = [oldf; newf]; ck = [oldc; newc];
    xk_extra = [oldxe_orig; newxe];
    % fprintf('distance fail new point is not adding \n');
          
else
    xk = x; fk = f; ck = c; xk_extra = x_extra; % nothing changes
end

end
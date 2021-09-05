
function newx  = boundary_check(newx, xu_bound, xl_bound, varargin)
if ~isempty(varargin)
    prob = varargin{1};
end


upper_check = xu_bound - newx;
if any(upper_check < 0) % beyond upper bound
    adjust_index = find(upper_check<0);
    newx(adjust_index) = xu_bound(adjust_index);
    
    if ~isempty(varargin)
        fprintf('problem %s upper bound violation \n', prob.name);
    end
end


lower_check = newx - xl_bound;
if any(lower_check < 0) % beyond lower bound
    adjust_index = find(lower_check<0);
    newx(adjust_index) = xl_bound(adjust_index);
    
    if ~isempty(varargin)
        fprintf('problem %s lower bound violation \n', prob.name);
    end
end   

end
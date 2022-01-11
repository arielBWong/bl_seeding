
function onboundary = onboundary_check(xl_startingpoint, xl_bl, xl_bu, bl, bu)
% when this function is called, we already know that xl is local optimal
onboundary = false;

xl_startingpoint = (xl_startingpoint - bl) ./ (bu - bl); 
xl_bl = (xl_bl - bl) ./ (bu - bl);
xl_bu = (xl_bu - bl) ./ (bu - bl);

if any( abs(xl_startingpoint - xl_bl) < 1e-5)
    onboundary = true;
    return;
end

if any(abs(xl_bu - xl_startingpoint) < 1e-5)
    onboundary = true;
    return
end

end

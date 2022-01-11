
function onboundary = onboundary_checkTrh(xl_startingpoint, xl_bl, xl_bu, trh, bl, bu)
% when this function is called, we already know that xl is local optimal
onboundary = false;

xl_startingpoint = (xl_startingpoint - bl) ./ (bu - bl); 
xl_bl = (xl_bl - bl) ./ (bu - bl);
xl_bu = (xl_bu - bl) ./ (bu - bl);


if any( abs(xl_startingpoint - xl_bl) < trh)
    onboundary = true;
    return;
end

if any(abs(xl_bu - xl_startingpoint) < trh)
    onboundary = true;
    return
end

end

function [xu, xl, fu, fl, fc, lc] = add_entry_extened(xu, xl, fu, fl, fc, lc, newxu, newxl, newfu, newfl, newfc, newlc)
xu  = [xu; newxu];
xl   = [xl; newxl];
fu  = [fu; newfu];
fl  = [fl; newfl];
fc  = [fc; newfc];
lc  = [lc; newlc];
end
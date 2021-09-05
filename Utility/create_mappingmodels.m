function [mapping_bl, mapping_upper, mapping_lower, mapping_lowvar]= create_mappingmodels(archive, param)



mapping_xu  = archive.xu(archive.flag == 1, :);
mapping_xl   = archive.xl(archive.flag == 1, :);
mapping_fu  = archive.fu(archive.flag == 1, :);
mapping_fl   = archive.fl(archive.flag == 1, :);


mapping_bl = Train_GPR(mapping_xu, mapping_fl, param);

mapping_x = [mapping_xu, mapping_xl];
mapping_upper = Train_GPR(mapping_x, mapping_fu, param);
mapping_lower = Train_GPR(mapping_x, mapping_fl, param);

mapping_lowvar    = Train_GPR(mapping_xu, mapping_xl, param);

end
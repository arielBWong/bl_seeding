function  [newxu, info] = Believer_nextExtended(archive,  prob, param_ea, krg_param, model_xl, arc_trgxl, extended, varargin)


%----------------------------------
%
if ~isempty(varargin)
    fighn    = varargin{1};
end

xu   = archive.xu;
xl    = archive.xl;
fu   = archive.fu;
fl     = archive.fl;
cu   = archive.cu;
cl    = archive.cl;

collectx = [xu,  xl];


% create model for infill (xu, xl) -> Fu
[~, of] = size(fu);
if of > 1
    fu_norm = normalization_nd(fu);
    ndf_indx = Paretoset(fu_norm);
    ndf = fu_norm(ndf_indx, :);
end

% ----------
% create infill search objective
if  of > 1
    if ~extended % original scheme
        arc_trg.x = xu;
        arc_trg.muf = fu_norm;
        
        krg_obj = Train_GPR(xu, fu_norm, krg_param);
        funh_obj  = @(xu)predict_kbhv(xu, krg_obj, ndf, krg_param, arc_trg);
        
        if ~isempty(cu)
            
            arc_trgc.x = xu;
            arc_trgc.muf =cu;
            funh_con       = @(xu)predict_cons(xu, krg_con, krg_param, arc_trgc);
        else
            funh_con  = @(xu)no_cons(xu);
        end
    else  % extended scheme
        
        arc_trg.x = [xu, xl];
        arc_trg.muf = fu_norm;
        
        krg_obj = Train_GPR([xu, xl], fu_norm, krg_param);
        funh_obj = @(xu)predict_kbhv(xu, krg_obj, ndf, krg_param, arc_trg, model_xl, arc_trgxl);
        if ~isempty(cu)
        
            arc_trgc.x = [xu, xl];
            arc_trgc.muf =cu;
            
            funh_con = @(xu)predict_cons(xu, krg_con, krg_param, arc_trgc, model_xl, arc_trgxl);
        else
            funh_con = @(xu)no_cons(xu);
        end
    end
else
    % single objective
    % objective
    if ~extended % original scheme
        % surrogate objectives
        arc_trg.x = xu;
        arc_trg.muf = fu;
        krg_obj = Train_GPR(xu,fu, krg_param);
        funh_obj  = @(xu)Predict_GPR(krg_obj, xu, krg_param, arc_trg);
        % surrogate constraints    
         if ~isempty(cu)
            krg_con = Train_GPR(xu, cu, krg_param);
            arc_trgc.x = xu;
            arc_trgc.muf =cu;
            funh_con       = @(xu)Predict_GPR(krg_con, xu, krg_param, arc_trgc);
        else
            funh_con  = @(xu)no_cons(xu);
        end
        
        
    else % extended
        x = [xu, xl];
        arc_trg.x = x;
        arc_trg.muf = fu;
        krg_obj  = Train_GPR(x,fu, krg_param);
        funh_obj  = @(xu)predict_extended(krg_obj, xu, krg_param, arc_trg, model_xl, arc_trgxl);
        
       
        if ~isempty(cu)
            arc_trgc.x = x;
            arc_trgc.muf = cu;
            krg_con = Train_GPR(x, cu, krg_param);
            funh_con = @(xu)predict_extendedcon(krg_con, xu, krg_param, arc_trgc, model_xl, arc_trgxl);
        else
            funh_con  = @(xu)no_cons(xu);
        end
        
    end
    
    
end

[~,~,~, archive]     = gsolver(funh_obj, prob.n_uvar, prob.xu_bl, prob.xu_bu, [], funh_con, param_ea);

newxu                    = archive.pop_last.X(1, :);
newxu                    = boundary_check(newxu, prob.xu_bu, prob.xu_bl);

info = [];

end

function [c] = no_cons(x)
c = [];
end

function f = predict_extended(krg_obj, xu, param, arc_trg, model_xl, arc_trgxl)
xl = Predict_GPR(model_xl, xu, param, arc_trgxl);
xx = [xu, xl];
f = Predict_GPR(krg_obj, xx, param, arc_trg);
end

function c = predict_extendedcon(krg_con, xu, param, arc_trg, model_xl, arc_trgxl)
xl = Predict_GPR(model_xl, xu, param, arc_trgxl);
xx = [xu, xl];
c = Predict_GPR(krg_con, xx, param, arc_trg);
end


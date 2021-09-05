function[best_x, info] = Believer_localsearch(krg_obj, krg_con, xu_bound, xl_bound, ...
                                 num_pop, num_gen, varargin)
   % fitness is an useless position, but cannot be deleted at this stage
                                                                               
    funh_obj = @(x)surrogate_predict(x, krg_obj, daceflag);
    funh_con = @(x)surrogate_predict(x, krg_con, daceflag);
    
    param.gen     = num_gen;
    param.popsize = num_pop;
    
    l_nvar           = size(xu_bound, 2);
    [~,~,~, archive] = gsolver(funh_obj, l_nvar,  xl_bound, xu_bound, [], funh_con, param);    
    best_x           = archive.pop_last.X(1, :);
    %-----stability
    info     = struct();
    info.krg = krg_obj;
    if  ~isempty(krg_con)
        info.krgc = krg_con;
    else
        info.krgc = [];
    end
    
    info.stable = true;
    
end




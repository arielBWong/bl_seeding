function [trgx_cheap, trgf_cheap] = nooverlap_fix(trgx_exp, trgf_exp, trgx_cheap, trgf_cheap)
% if expensive solutions are not included in cheap solutions
% augment cheap solutions
lia = ismember(trgx_exp, trgx_cheap);
if sum(lia) == 0 % means no overlap
    k = oodacefit(trgx_cheap,trgf_cheap); % if input is not cell means kriging
    
    % predict expensive in cheap f
    [exp_prd,~] = k.predict(trgx_exp);
    
    % augment x
    trgx_cheap = [trgx_cheap; trgx_exp];
    trgf_cheap = [trgf_cheap; exp_prd];   
end

end
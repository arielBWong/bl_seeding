function[co_mdl, co_mdlparam] = cokrgmodel_building(x_trg, y_trg, prob)
% this function build cokrg model with f being normalised

% retrieve training data from input
            
cheap_f = y_trg{1,:};
nf_cheap = size(cheap_f, 1);

expensive_f = y_trg{2, :};
nf_expensive = size(expensive_f, 1);

collection_f = [cheap_f; expensive_f];

[collection_fnorm, mu, sigma] = normalization_z(collection_f);
co_mdlparam.f_mu = mu;
co_mdlparam.f_sigma = sigma;

y_trg{1, :} = collection_fnorm(1:nf_cheap, :);
y_trg{2, :} = collection_fnorm(nf_cheap+1:end, :);


x_trg{1, :} = (x_trg{1, :} - prob.xl_bl)./ (prob.xl_bu - prob.xl_bl);
x_trg{2, :} = (x_trg{2, :} - prob.xl_bl)./ (prob.xl_bu - prob.xl_bl);

co_mdl = oodacefit(x_trg, y_trg);

end


function [f_norm, mu, sigma] = normalization_z(f)
% single normalization with zscore
if size(f, 2)> 1
    error('zscore is not for mo');
end
[f_norm, mu, sigma] = zscore(f, 0, 1); 
end


function f = denormzscore(f, mu, sigma)

f = f * sigma + mu;
end

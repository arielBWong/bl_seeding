function cov = covariation(x, mu)
% https://au.mathworks.com/help/matlab/ref/cov.html
% calculate covariance for multivariant normal distribution assumption

[n, m] = size(x);
cov = zeros(m, m);
for i = 1:m
    for j = i:m
        
        two_columns = [x(:, i), x(:, j)];
        two_means = [mu(i), mu(j)];
        cov_ij = sum(prod( (two_columns - two_means), 2))/(n-1);
        cov(i, j) = cov_ij;
    end 
end

var_in_diag = diag(diag(cov)); % create matrix with only diagnal of cov
cov_upperTri = cov - var_in_diag;
cov_lowerTri = transpose(cov_upperTri);

cov = cov_upperTri + var_in_diag + cov_lowerTri;

end
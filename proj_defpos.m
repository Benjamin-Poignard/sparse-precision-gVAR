function M_pos = proj_defpos(M)

% the following step is proposed page 2403 in Cai and Zhou (2012),
% "Optimal rates of convergence for sparse covariance matrix estimation",
% AoS, Vol. 40, No. 5, 2389-2420
p = size(M,2);
[V,D] = eig(M); eigenvalues = diag(D);
eigenvalues(eigenvalues<0)=0; M_pos = zeros(p);
for i = 1:p
    M_pos = M_pos + eigenvalues(i)*V(:,i)*V(:,i)';
end

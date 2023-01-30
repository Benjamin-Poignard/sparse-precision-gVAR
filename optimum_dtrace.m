function Sol = optimum_dtrace(A,B)

% Inputs: - A: symmetric positive-definite matrix
%         - B: any q-dimensional matrix

% Output: -Sol: solution provided in Theorem 1 in Zhang and Zou (2014), 
%          'Sparse precision matrix estimation via lasso penalized D-trace 
%          loss', Biometrika 

[V,D] = eigs(A,size(A,1)); C = zeros(size(A,2));
for i=1:size(A,2)
   for j = i:size(A,2)
       C(i,j) = 2/(D(i,i)+D(j,j));
       C(j,i) = C(i,j);
   end
end
Btemp = V'*B*V;
Sol = V*(Btemp.*C)*V';

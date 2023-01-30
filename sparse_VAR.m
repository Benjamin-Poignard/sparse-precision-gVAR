function b = sparse_VAR(data,p,lambda,method)

% Inputs: - X: T x N matrix of data
%         - p: number of lags in the VAR process
%         - lambda: regularization parameter (grid)
%         - method: 'lasso', 'ridge'

% Output: - b: estimated sparse VAR parameter (autoregressive p x p matrix
%           paramter

[T,~] = size(data);

% creation of the vector of covariate
X = []; data = data';
for tt = p+1:T
    x_temp_reg = [];
    for kk = 1:p
        x_temp_reg = [x_temp_reg ; data(:,tt-kk)];
    end
    X = [X , x_temp_reg];
end
X = X'; data = data';
b = [];
for ii = 1:size(data,2)
    [b1,~] = VAR_penalized(data(p+1:end,ii),X,lambda,method);
    b = [b;b1'];
end
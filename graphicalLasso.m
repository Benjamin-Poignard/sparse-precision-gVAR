function [theta,W] = graphicalLasso(X,rho,maxIt,tol)

% Graphical Lasso procedure proposed by Friedman et al. (2007), 'Sparse 
% inverse covariance estimation with the graphical lasso', Biostatistics.
% Minimize the loss: tr(S*Theta) - log|Theta| + rho * ||Theta||_1

% Note: This function uses the function *lassoShooting* (shooting
% algorithm) to solve the LASSO problem. But any alternative LASSO 
% algorithm may also be considered.

% Input: - X: sample of observations
%        - rho: tuning parameter
%        - maxIt: maximum number of iterations (optional)
%        - tol: convergence tolerance level (optional)

% Outputs: - theta: vech(Theta), i.e. column vector that stacks the columns
%            of the lower triangular part of Theta with Theta the sparse
%            precision matrix estimated by the LASSO Gaussian loss
%          - W: regularized covariance matrix estimator, W = Theta^-1

S = cov(X); p = size(S,1);

if nargin < 4, tol = 1e-6; end
if nargin < 3, maxIt = 1e5; end

% Initialization
W = S;   % diagonal of W remains unchanged, no penalization for diagonal elements
W_old = W;
i = 0;

% Graphical Lasso loop
while i < maxIt
    i = i+1;
    for j = p:-1:1
        jminus = setdiff(1:p,j);
        [V,D] = eig(W(jminus,jminus));
        d = diag(D);
        X = V * diag(sqrt(d)) * V'; % W_11^(1/2)
        Y = V * diag(1./sqrt(d)) * V' * S(jminus,j);    % W_11^(-1/2) * s_12
        b = lassoShooting(X,Y,rho,maxIt,tol);
        W(jminus,j) = W(jminus,jminus) * b;
        W(j,jminus) = W(jminus,j)';
    end
    % Stop criterion
    if norm(W-W_old,1) < tol
        break;
    end
    W_old = W;
end
if i == maxIt
    fprintf('%s\n', 'Maximum number of iteration reached, glasso may not converge.');
end

Theta = W^-1; Theta(abs(Theta)<0.0001)=0; theta = vech(Theta);

% Shooting algorithm for Lasso (unstandardized version)
function b = lassoShooting(X,Y,lambda,maxIt,tol)

if nargin < 5, tol = 1e-6; end
if nargin < 4, maxIt = 1e2; end

% Initialization
[n,p] = size(X);
if p > n
    b = zeros(p,1); % From the null model, if p > n
else
    b = X \ Y;  % From the OLS estimate, if p <= n
end
b_old = b;
i = 0;

% Precompute X'X and X'Y
XTX = X'*X;
XTY = X'*Y;

% Shooting loop
while i < maxIt
    i = i+1;
    for j = 1:p
        jminus = setdiff(1:p,j);
        S0 = XTX(j,jminus)*b(jminus) - XTY(j);  % S0 = X(:,j)'*(X(:,jminus)*b(jminus)-Y)
        if S0 > lambda
            b(j) = (lambda-S0) / norm(X(:,j),2)^2;
        elseif S0 < -lambda
            b(j) = -(lambda+S0) / norm(X(:,j),2)^2;
        else
            b(j) = 0;
        end
        
    end
    delta = norm(b-b_old,1);    % Norm change during successive iterations
    if delta < tol, break; end
    b_old = b;
end
if i == maxIt
    fprintf('%s\n', 'Maximum number of iteration reached, shooting may not converge.');
end
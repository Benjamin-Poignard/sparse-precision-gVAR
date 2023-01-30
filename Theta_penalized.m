function theta = Theta_penalized(X,loss,lambda,method,a_scad,b_mcp)

% Inputs: - X: vector of observation
%         - loss: 'Gaussian' or 'DTrace'; if 'Gaussian', the GLASSO
%         algorithm is performed (only with LASSO)
%         - lambda: regularization parameter
%         - method: 'scad', 'mcp', 'lasso'
%         - a_scad: value of the scad parameter
%         - b_mcp: value of the mcp parameter

% Output: - theta: vech(Theta), i.e. column vector that stacks the
%           columns of the lower triangular part of Theta

switch loss
    case 'DTrace'
        theta = sparse_dtrace(X,lambda,method,a_scad,b_mcp);
    case 'Gaussian'
        [theta,~] = graphicalLasso(X,lambda); 
end
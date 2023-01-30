function theta = sparse_dtrace(X,lambda,method,a_scad,b_mcp)

% Input: - X: sample of observations
%        - lambda: tuning parameter
%        - method: 'lasso', 'scad', 'mcp'
%        - a_scad: value of the scad parameter
%        - b_mcp: value of the mcp parameter

% Output: - theta: vech(Theta), i.e. column vector that stacks the columns
%           of the lower triangular part of Theta

% For the SCAD and MCP penalty functions, the local linear approximation
% (LLA) algorithm of Zou and Li (2008), 'One-step sparse estimates in
% non-concave penalized likelihood models', the Annals of Statistics, is
% implemented.
% Note: rather than updating the weights as in problem (2.7) of Zou and Li
% (2008), the standard LASSO is run first; then the entries of this
% estimator enter the weight of (2.7) in Zou and Li (2008), where the
% weight is the first order derivative of the SCAD/MCP evaluated at the
% LASSO estimator.

p = size(X,2); S = cov(X);
if min(eig(S))<eps
    S = proj_defpos(S);
end

switch method
    
    case 'lasso'
        
        Theta = dtrace_admm(S,eye(p),lambda,method,a_scad,b_mcp); theta = vech(Theta);
        
    case 'alasso'
        
        if p<size(X,1)
            Theta = dtrace_admm(S,inv(S),lambda,method,a_scad,b_mcp); theta = vech(Theta);
        else
            Theta_lasso = dtrace_admm(S,eye(p),lambda,'lasso',a_scad,b_mcp);
            Theta = dtrace_admm(S,Theta_lasso,lambda,method,a_scad,b_mcp); theta = vech(Theta);
        end
        
        
    case 'scad'
        
        % first step: standard LASSO estimation
        Theta_lasso = dtrace_admm(S,eye(p),lambda,'lasso',a_scad,b_mcp);
        % second step: use Theta_lasso to perform the LLA algorithm
        Theta = dtrace_admm(S,Theta_lasso,lambda,method,a_scad,b_mcp); theta = vech(Theta);
        
    case 'mcp'
        
        % first step: standard LASSO estimation
        Theta_lasso = dtrace_admm(S,eye(p),lambda,'lasso',a_scad,b_mcp);
        % second step: use Theta_lasso to perform the LLA algorithm
        Theta = dtrace_admm(S,Theta_lasso,lambda,method,a_scad,b_mcp); theta = vech(Theta);
        
end
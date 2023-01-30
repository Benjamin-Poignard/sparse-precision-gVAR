function [theta_est,lambda_opt]= sparse_precision(X,loss,lambda,method,a_scad,b_mcp)

% Inputs: - X: vector of observation
%         - loss: 'Gaussian' or 'DTrace'; if 'Gaussian', the GLASSO
%         algorithm is performed (only with LASSO)
%         - lambda: tuning parameter
%         - method (for D-trace loss): 'scad', 'mcp', 'alasso', 'lasso'
%         - a_scad: value of the scad parameter
%         - b_mcp: value of the mcp parameter

% Outputs: - theta_est: vech(Theta), i.e. column vector that stacks the
%           columns of the lower triangular part of Theta
%          - lambda_opt: optimal tuning parameter value chosen by
%          cross-validation; if no-cross validation is performed (i.e., the
%          input lambda is a scalar), then lambda_opt = lambda

[T,d] = size(X); dim = d*(d+1)/2;

% if lambda is a vector of tuning parameters, then apply cross-validation
if length(lambda)>1
    % construction of the train sets and validation datasets
    len_in = round(0.75*T); X_in = X(1:len_in,:); theta_fold = zeros(dim,length(lambda));
    X_out = X(len_in+1:T,:);
    parfor jj = 1:length(lambda)
        theta_fold(:,jj) = Theta_penalized(X_in,loss,lambda(jj),method,a_scad,b_mcp);
    end
    % evaluate the loss over the test datasets using the estimators
    % obtained over the train datasets
    count = zeros(length(lambda),1);
    for ii = 1:length(lambda)
        Theta = dvech(theta_fold(:,ii),d); S = cov(X_out);
        switch loss
            case 'Gaussian'
                L = trace(S*Theta)-log(det(Theta));
            case 'DTrace'
                L =  0.5*trace(Theta^2*S)-trace(Theta);
        end
        count(ii) = count(ii) + L;
    end
    clear ii kk
    
    ii = count==min(min(count)); lambda_opt = lambda(ii);
    if length(lambda_opt)>1
        lambda_opt = lambda(1);
    end
    theta_est = Theta_penalized(X,loss,lambda_opt,method,a_scad,b_mcp);
else
    lambda_opt = lambda;
    % if lambda is a scalar, no cross-validation is performed
    theta_est = Theta_penalized(X,loss,lambda_opt,method,a_scad,b_mcp);
end
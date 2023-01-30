function [theta_est,lambda_opt]= VAR_penalized(Y,X,lambda,method)

% Inputs : - Y: response variable
%          - X: covariates
%          - lambda: tuning parameter
%          - method: 'ridge' or 'lasso'

% Outputs: - theta_est: penalized solution
%          - lambda_opt: optimal tuning parameter value chosen by
%           cross-validation; if no-cross validation is performed (i.e. the
%           input lambda is a scalar), then lambda_opt = lambda

[T,d] = size(X);

% if lambda is a vector of tuning parameters, then apply cross-validation
if length(lambda)>1
    % construction of the train sets and validation datasets
    len_in = round(0.75*T);
    X_in = X(1:len_in,:); X_out = X(len_in+1:T,:);
    Y_in = Y(1:len_in,:); Y_out = Y(len_in+1:T,:);
    theta_fold = zeros(d,length(lambda));
    
    parfor jj = 1:length(lambda)
        switch method
            case 'lasso'
                theta_fold(:,jj) = lasso(X_in,Y_in,'Lambda',lambda(jj),'Intercept',false);
            case 'ridge'
                theta_fold(:,jj) = ridge(Y_in,X_in,lambda(jj));
        end
    end
    
    % evaluate the loss over the test dataset using the estimators
    % obtained over the train dataset
    count = zeros(length(lambda),1);
    for ii = 1:length(lambda)
        L = sum((Y_out-X_out*theta_fold(:,ii)).^2);
        count(ii) = count(ii) + L;
    end
    clear ii
    
    ii = count==min(min(count)); lambda_opt = lambda(ii);
    if length(lambda_opt)>1
        lambda_opt = lambda(1);
    end
    switch method
        case 'lasso'
            theta_est = lasso(X,Y,'Lambda',lambda_opt,'Intercept',false);
        case 'ridge'
            theta_est = ridge(Y,X,lambda_opt);
    end
else
    lambda_opt = lambda;
    % if lambda is a scalar, no cross-validation is performed
    switch method
        case 'lasso'
            theta_est = lasso(X,Y,'Lambda',lambda_opt,'Intercept',false);
        case 'ridge'
            theta_est = ridge(Y,X,lambda_opt);
    end
end
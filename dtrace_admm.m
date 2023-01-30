function Theta = dtrace_admm(S,W,lambda,method,a_scad,b_mcp)

% ADMM algorithm proposed in Zhang and Zou (2014), 'Sparse precision matrix
% estimation via lasso penalized D-trace loss', Biometrika
% The main code corresponds to Algorithm 2 of Zhang and Zou (2014); if the
% estimated Theta is not positive definite, then Algorithm 1 of Zhang and
% Zou (2014) is performed.

% Inputs: - S: sample variance covariance matrix of the observations
%         - W: matrix of weights active for the local linear approximation
%         algorithm for SCAD and MCP; for SCAD and MCP, W is the LASSO
%         estimator; for the LASSO penalization, W = eye(p)
%         - lambda: tuning parameter
%         - method: 'scad', 'mcp', 'alasso', 'lasso'
%         - a_scad: value of the scad parameter
%         - b_mcp: value of the mcp parameter

% Output: - Theta: sparse precision matrix estimated by the D-trace loss

% Algorithm 2 of Zhang and Zou (2014)
p = size(S,2); Theta0(:,:,1) = inv(diag(diag(S))); Lambda(:,:,1) = Theta0(:,:,1);
Theta(:,:,1)=zeros(p);
maxIt = 1e5; crit = eps; k = 1; i = 0;
while i < maxIt
    i=i+1;
    k=k+1;
    Theta(:,:,k) = optimum_dtrace(S+eye(p),eye(p)+Theta0(:,:,k-1)-Lambda(:,:,k-1));
    Theta0(:,:,k) = soft_threshold_precision(Theta(:,:,k)+Lambda(:,:,k-1),W,lambda,method,a_scad,b_mcp);
    Lambda(:,:,k) = Lambda(:,:,k-1)+Theta(:,:,k)-Theta0(:,:,k);
    error1 = (norm(Theta(:,:,k)-Theta(:,:,k-1))^2/max([1,norm(Theta(:,:,k)),norm(Theta(:,:,k-1))]));
    error2 = (norm(Theta0(:,:,k)-Theta0(:,:,k-1))^2/max([1,norm(Theta0(:,:,k)),norm(Theta0(:,:,k-1))]));
    if ((error1<crit)&&(error2<crit))
        break;
    end
    
end
if i == maxIt
    fprintf('%s\n', 'Maximum number of iterations reached, ADMM may not converge.');
end
Theta = Theta(:,:,end);
% Verify the positive-definiteness of the solution Theta
% If not positive-definite, Algorithm 1 of Zhang and Zou (2014)
if min(eig(Theta))<0.0001
    Theta0(:,:,1) = Theta; Theta1(:,:,1) = Theta; Lambda0(:,:,1) = Theta0(:,:,1); Lambda1(:,:,1) = Theta1(:,:,1);
    crit = eps; k = 1; i = 0;
    while i < maxIt
        i=1+1;
        k=k+1;
        Theta(:,:,k) = optimum_dtrace(S+2*eye(p),eye(p)+Theta0(:,:,k-1)-Lambda0(:,:,k-1)-Lambda1(:,:,k-1));
        Theta0(:,:,k) = soft_threshold_precision(Theta(:,:,k)+Lambda0(:,:,k-1),W,lambda,method,a_scad,b_mcp);
        Theta1(:,:,k) = proj_defpos(Theta(:,:,k)+Lambda1(:,:,k-1));
        Lambda0(:,:,k) = Lambda0(:,:,k-1)+Theta(:,:,k)-Theta0(:,:,k);
        Lambda1(:,:,k) = Lambda1(:,:,k-1)+Theta(:,:,k)-Theta1(:,:,k);
        error1 = (norm(Theta(:,:,k)-Theta(:,:,k-1))^2/max([1,norm(Theta(:,:,k)),norm(Theta(:,:,k-1))]));
        error2 = (norm(Theta0(:,:,k)-Theta0(:,:,k-1))^2/max([1,norm(Theta0(:,:,k)),norm(Theta0(:,:,k-1))]));
        if ((error1<crit)&&(error2<crit))
            break;
        end
    end
end
if i == maxIt
    fprintf('%s\n', 'Maximum number of iterations reached, ADMM may not converge.');
end
Theta(abs(Theta)<10^(-7))=0;
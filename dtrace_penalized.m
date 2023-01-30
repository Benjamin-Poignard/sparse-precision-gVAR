function theta = dtrace_penalized(X,lambda,R,method,a_scad,b_mcp)

p = size(X,2); dim = p*(p+1)/2;
switch method
    case 'scad'
        mu = 1/(a_scad-1);
    case 'mcp'
        mu = 1/b_mcp;
    case 'lasso'
        mu = 0;
end
eta = 5000;
% calibration of the step size
nu = (1/eta)/(1+(mu/eta));

optimoptions.MaxRLPIter = 200000;
optimoptions.MaxFunEvals = 200000;
optimoptions.Algorithm = 'sqp';
optimoptions.TolCon = 1e-09;
optimoptions.TolRLPFun = 1e-09;
optimoptions.MaxSQPIter = 200000;
optimoptions.Display = 'iter';

max_iter = 100;
param = zeros(dim,1);

S = cov(X);
if min(eig(S))<eps
    S = proj_defpos(S);
end
Sinv = inv(S); param_update = vech(Sinv);

count = 0;
while (norm(param_update - param)^2/max([1,norm(param_update),norm(param)]))>eps
    
    param = param_update;
    Theta = dvech(param,p);
    
    gradient = vech((Theta*S+S*Theta)/2-eye(p));
    gradient_modified = gradient - mu*param;
    Z = (1/(1+(mu/eta)))*(param-gradient_modified/eta);
    param_est = zeros(dim,1);
    switch method
        case 'scad'
            cnt = 1;
            for ii = 1:dim
                tmplambda = lambda;
                if (ii == (cnt-1)*(p+1-cnt/2)+1)
                    tmplambda = 0;
                    cnt = cnt + 1;
                end
                if (0 <= abs(Z(ii)) && abs(Z(ii)) <= nu*tmplambda)
                    param_est(ii) = 0;
                elseif (nu*tmplambda <= abs(Z(ii)) && abs(Z(ii)) <= (nu+1)*tmplambda)
                    param_est(ii) = Z(ii)-(sign(Z(ii))*nu*tmplambda);
                elseif ((nu+1)*tmplambda <= abs(Z(ii)) && abs(Z(ii)) <= a_scad*tmplambda)
                    param_est(ii) = (Z(ii)-((sign(Z(ii))*a_scad*nu*tmplambda)/(a_scad-1)))/(1-(nu/(a_scad-1)));
                elseif (a_scad*tmplambda <= abs(Z(ii)))
                    param_est(ii) = Z(ii);
                end
            end
            if (side_constraint(param_est,lambda,method,a_scad,b_mcp) <= R)
                param_update = param_est;
            else
                count = count + 1;
                if (count > max_iter)
                    break
                end
                start = vech(diag(diag(Sinv)));
                [param_est,~,~,~,~,~]=fmincon(@(x)OLS_loss(x,param-gradient_modified/eta),vech(start),[],[],[],[],[],[],@(x)constr_R(x,p,method,lambda,a_scad,b_mcp,R),optimoptions);
                param_update = param_est; param_update(abs(param_update)<0.001)=0;
            end
        case 'mcp'
            cnt = 1;
            for ii = 1:dim
                tmplambda = lambda;
                if (ii == (cnt-1)*(p+1-cnt/2)+1)
                    tmplambda = 0;
                    cnt = cnt + 1;
                end
                if (0 <= abs(Z(ii)) && abs(Z(ii)) <= nu*tmplambda)
                    param_est(ii) = 0;
                elseif (nu*tmplambda <= abs(Z(ii)) && abs(Z(ii)) <= b_mcp*tmplambda)
                    param_est(ii) = (Z(ii)-(sign(Z(ii))*nu*tmplambda))/(1-nu/b_mcp);
                elseif (b_mcp*tmplambda <= abs(Z(ii)))
                    param_est(ii) = Z(ii);
                end
            end
            if (side_constraint(param_est,lambda,method,a_scad,b_mcp) <= R)
                param_update = param_est;
            else
                count = count + 1;
                if (count > max_iter)
                    break
                end
                start = vech(diag(diag(Sinv)));
                [param_est,~,~,~,~,~]=fmincon(@(x)OLS_loss(x,param-gradient_modified/eta),vech(start),[],[],[],[],[],[],@(x)constr_R(x,p,method,lambda,a_scad,b_mcp,R),optimoptions);
                param_update = param_est; param_update(abs(param_update)<0.001)=0;
            end
        case 'lasso'
            cnt = 1;
            for ii = 1:dim
                tmplambda = lambda;
                if (ii == (cnt-1)*(p+1-cnt/2)+1)
                    tmplambda = 0;
                    cnt = cnt + 1;
                end
                param_est(ii) = sign(Z(ii))*subplus(abs(Z(ii))-tmplambda/eta);
            end
            if (side_constraint(param_est,lambda,method,a_scad,b_mcp) <= R)
                param_update = param_est;
            else
                count = count + 1;
                if (count > max_iter)
                    break
                end
                start = vech(diag(diag(Sinv)));
                [param_est,~,~,~,~,~]=fmincon(@(x)OLS_loss(x,param-gradient_modified/eta),vech(start),[],[],[],[],[],[],@(x)constr_R(x,p,method,lambda,a_scad,b_mcp,R),optimoptions);
                param_update = param_est; param_update(abs(param_update)<0.001)=0;
            end
    end
    count = count+1;
end
theta = param_update;

% final step: verifying the positive-definiteness
if min(eig(dvech(theta,p)))<eps
    Theta_pos = proj_defpos(dvech(theta,p));
    theta = vech(Theta_pos);
end
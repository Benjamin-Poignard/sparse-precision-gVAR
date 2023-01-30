function Solpen = soft_threshold_precision(A,Theta_weight,lambda,method,a_scad,b_mcp)

% Inputs: - A: symmetric positive-definite matrix to be penalized
%         - Theta_weight: initial solution for LLA algorithm defined as the 
%         LASSO for SCAD and MCP and the inverse of the sample covariance
%         matrix for the adaptive LASSO
%         - lambda: regularization parameter
%         - method (for D-trace loss): 'scad', 'mcp', 'alasso', 'lasso'
%         - a_scad: value of the scad parameter
%         - b_mcp: value of the mcp parameter

% Output: - Solpen: penalized solution

p = size(A,2); dim = p*(p+1)/2; Z = vech(A); Solpen = zeros(dim,1); cnt = 1; t_weight = vech(Theta_weight);
for i = 1:dim
    tmplambda = lambda;
    if (i == (cnt-1)*(p+1-cnt/2)+1)
        tmplambda = 0;
        cnt = cnt + 1;
    end
    if (tmplambda == 0)
        Solpen(i) = Z(i);
    else
        switch method
            case 'lasso'
                if (Z(i)>tmplambda)
                    Solpen(i) = Z(i)-tmplambda;
                elseif (Z(i)<-tmplambda)
                    Solpen(i) = Z(i)+tmplambda;
                elseif (-tmplambda <= Z(i) && Z(i) <= tmplambda)
                    Solpen(i) = 0;
                end
            case 'alasso'
                gamma = 0.5;
                if (Z(i)>tmplambda*(abs(t_weight(i))^(-gamma)))
                    Solpen(i) = Z(i)-tmplambda*(abs(t_weight(i))^(-gamma));
                elseif (Z(i)<-tmplambda*(abs(t_weight(i))^(-gamma)))
                    Solpen(i) = Z(i)+tmplambda*(abs(t_weight(i))^(-gamma));
                elseif (-tmplambda*(abs(t_weight(i))^(-gamma)) <= Z(i) && Z(i) <= tmplambda*(abs(t_weight(i))^(-gamma)))
                    Solpen(i) = 0;
                end
            case 'scad'
                if (abs(t_weight(i))<=tmplambda)
                    scad_dev = tmplambda;
                elseif (tmplambda<abs(t_weight(i)))
                    scad_dev = subplus(a_scad*tmplambda-abs(t_weight(i)))/(a_scad-1);
                end
                if (Z(i)>scad_dev)
                    Solpen(i) = Z(i)-scad_dev;
                elseif (Z(i)<-scad_dev)
                    Solpen(i) = Z(i)+scad_dev;
                elseif (-scad_dev <= Z(i) && Z(i) <= scad_dev)
                    Solpen(i) = 0;
                end
            case 'mcp'
                if (abs(t_weight(i))<=b_mcp*tmplambda)
                    mcp_dev = (b_mcp*tmplambda-abs(t_weight(i)))/b_mcp;
                elseif (abs(t_weight(i))>b_mcp*tmplambda)
                    mcp_dev = 0;
                end
                if (Z(i)>mcp_dev)
                    Solpen(i) = Z(i)-mcp_dev;
                elseif (Z(i)<-mcp_dev)
                    Solpen(i) = Z(i)+mcp_dev;
                elseif (-mcp_dev <= Z(i) && Z(i) <= mcp_dev)
                    Solpen(i) = 0;
                end
        end
    end
end
Solpen = dvech(Solpen,p);
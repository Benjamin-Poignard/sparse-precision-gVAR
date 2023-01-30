function N = check_NZ(beta0,beta)

% Input:
%       - beta_0: true vector of parameters 
%       - beta: estimated vector of parameters
% Output:
%       - N: number of correctly estimated non-zero coefficients 


p = length(beta); N=0;
for i = 1:p
   if (abs(beta(i))>0)&&(abs(beta0(i))>0)
	N=N+1;
   end
end
function C = check_zero(beta0,beta)

% Input:
%       - beta_0: true vector of parameters 
%       - beta: estimated vector of parameters
% Output:
%       - N: number of correctly estimated zero coefficients 

p = length(beta);
C=0;
for i = 1:p
   if beta0(i)==0
       if beta(i)==beta0(i)
           C=C+1;
       end
   end
end
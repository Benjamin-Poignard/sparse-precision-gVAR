function [N1,N2] = check_IC_2(beta0,beta)

% Input:
%       - beta_0: true vector of parameters 
%       - beta: estimated vector of parameters
% Output:
%       - N1: number of true non-zero coefficients incorrectly estimated, i.e.
%       estimated zero coefficients whereas the true coefficient is a non-zero
%       - N2: number of true zero coefficients incorrectly estimated, i.e.
%        estimated non-zero coefficients whereas the true coefficient is a zero 

p = length(beta);
N1=0; N2=0;
for i = 1:p
   if (beta(i)==0)&&(abs(beta0(i))>0)
	N1=N1+1;
   end
   if (abs(beta(i))>0)&&(beta0(i)==0)
        N2=N2+1;
   end
end
function  [iMS,iRS,iMAE] = sparsistency_check_general(vP_true,vP_est)
  % Check sparsistency and consistency
  iP = length(vP_true);
  vC = (abs(vP_true) > 0); vC1 = (1:iP)'.*vC; iPC = sum(vC);
  vS = 1 - vC; vS1=(1:iP)'.*vS; iPS = sum(vS);
  vC0 = zeros(iPC,1); vS0 = zeros(iPS,1);j=0;h=0;
  for i=1:iP
      if (vC1(i)>0)
          j=j+1;
          vC0(j) = vC1(i);
      else
          h=h+1;
          vS0(h) = vS1(i);
      end
  end 
  
  iNS = length(vS0);
  %iMS = iNS - sum((abs(vP_est(vS0))>0));
  iMS = sum(abs(vP_est(vS0))<1e-06);
  iRS = iMS/iNS;
  iMAE = mean(abs((vP_true(vC0)-vP_est(vC0))));
  %[iMS iRS iMAE]
end
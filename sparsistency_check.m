function  [iMS,iRS,iMAE] = sparsistency_check(Theta_true,Theta_est,vS,vC)
  % Check sparsistency and consistency
  %
  % vech(Theta_true)
  % vS shows the locations such that the elements of vech(Theta_true) are
  % zeros.
  % vC shows the locations such that the elements of vech(Theta_true) are
  % non-zeros.
 
  vTheta_true = vech(Theta_true);
  vTheta_est = vech(Theta_est);
  
  iNS = length(vS);
  iMS = iNS - sum((abs(vTheta_est(vS))>0));
  iRS = iMS/iNS;
  iMAE = mean(abs((vTheta_true(vC)-vTheta_est(vC))));
  %[iMS iRS iMAE]
end
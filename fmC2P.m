function mPR = fmC2P(mC)
  % mC: covariance
  % mP: Partial autocorrelations
  mR = corrcov(mC);%mR
  iN = length(mR);
  imR = inv(mR);%imR
  mD = inv(diag(sqrt(diag(imR))));  
  mP = -mD*imR*mD;
  mPR = mP - diag(diag(mP)) + eye(iN);
  %vPS = ones(iN,1)./sqrt(diag(inv(mC)));
end
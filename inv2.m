function imS = inv2(mS)
  % inverse of symmetric matrix
  mL = chol(mS,'lower');
  imL = inv(mL);
  imS = imL'*imL;
end
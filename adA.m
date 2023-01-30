function mAp = adA(mA)
  [iN,iNP]=size(mA);
  iP = iNP/iN;
  mAp = zeros(iN,iN*iP);
  for i=1:iP
      mAp(:,(i-1)*iN+1:i*iN) = mA(:,(iP-i)*iN+1:(iP+1-i)*iN);
  end
end
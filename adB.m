function mB = adB(mBp)
  [iN,iNP1]=size(mBp);
  iP = (iNP1-iN)/iN;
  mB = zeros(iN,iN*(iP+1));
  for i=1:(iP+1)
      mB(:,(i-1)*iN+1:i*iN) = mBp(:,(iP+1-i)*iN+1:(iP+2-i)*iN);
  end
end
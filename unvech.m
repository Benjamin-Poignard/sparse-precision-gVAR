function mX = unvech(vX)
%need to check
  iN = round(0.5*(-1+sqrt(1+8*length(vX))));%[length(vX) iN]
  mX = zeros(iN,iN);
  mX(:,1) = vX(1:iN);
  idx=iN+1;
  for j=2:iN;
      mX(j:end,j) = vX(idx:idx+(iN-j));
      idx = idx + iN-j+1;
  end
  mXu = triu(mX',1);
  mX = mX + mXu;
  
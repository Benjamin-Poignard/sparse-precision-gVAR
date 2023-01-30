clear
iN=20; iP=3;

rng(123,'twister');

% set parameters, mBp = [ B_0, B_1, ..., B_p]
mBp = SetB(iN,iP);

% set parameters for generate y_t
mB0 = mBp(:,1:iN);imB0 = inv(mB0);
mA = imB0*mBp(:,iN+1:end);mS_u = imB0*imB0';

% change ordering as mB = [ B_p, ..., B_1,B_0]
mB = adB(mBp);

% convert B_0,B_1,...B_p to \Sigma_x
mSig_x = fVAR2Sig(mB);
% Note: inv(mSig_x) is sparse

% equation (3) and (4)
imL = chol(mSig_x,'lower');
mL = inv(imL);
mB0s = mL(iN*iP+1:end,iN*iP+1:end);
mBas = -mL(iN*iP+1:end,1:iN*iP);
mBs = [mBas mB0s];

% check the difference
vD = reshape(mBs - mB,iN*iN*(iP+1),1);
iD = vD'*vD;iD

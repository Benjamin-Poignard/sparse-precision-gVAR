function [mBh,mB0h,mAh,mSu]=Theta2Ba(mTheta,iN)

% Inputs: - mTheta: precision matrix for the gVAR process
%         - iN: dimension of the random vector

% Outputs: - mBh: matrix coefficients of the SVAR (B1, B2, B3, ...)
%          - mB0h: matrix B0 of the SVAR
%          - mAh: matrix coefficients of the VAR (A1, A2, ...)
%          - mSu: variance-covariance matrix of the VAR system innovation


mSx = inv(mTheta);
d = length(mTheta); iP = length(mTheta)/iN-1;
%'Estimate B_p,...,B_1'
mPR_x = fmC2P(mSx);
vS_x = sqrt(diag(mSx));
mAd = zeros(d,d);
for j=1:iN
    r= d-iN+j;
for i=1:d
    mQ = mSx;
    mQ(i,:)=[];mQ(:,i)=[];
    if (i<r)
        mXX = mQ; mXX(r-1,:)=[];mXX(:,r-1)=[];
        vXY = mQ(r-1,:)'; vXY(r-1)=[];
        beta = mXX\vXY;
        iR2 = (beta'*mXX*beta)/mSx(r,r);
    elseif (i>r)
        mXX = mQ; mXX(r,:)=[];mXX(:,r)=[];
        vXY = mQ(r,:)'; vXY(r)=[];
        beta = mXX\vXY;
        iR2 = (beta'*mXX*beta)/mSx(r,r);
    else
        iR2=0;
    end
    mQ = mSx;
    mQ(r,:)=[];mQ(:,r)=[];
    if (i<r)
        mXX = mQ; mXX(i,:)=[];mXX(:,i)=[];
        vXY = mQ(i,:)'; vXY(i)=[];
        beta = mXX\vXY;
        iR2d = (beta'*mXX*beta)/mSx(i,i);
    elseif (i>r)
        mXX = mQ; mXX(i-1,:)=[];mXX(:,i-1)=[];
        vXY = mQ(i-1,:)'; vXY(i-1)=[];
        beta = mXX\vXY;
        iR2d = (beta'*mXX*beta)/mSx(i,i);
    else
        iR2d=0;
    end
    mAd(r,i)=sqrt((1-iR2)/(1-iR2d));        
end
vS_xd = kron(ones(iP+1,1),vS_x(end-iN+1:end));
mC_x = diag(vS_xd)*(mPR_x.*mAd)*inv(diag(vS_x));
mBsh = mC_x(end-iN+1:end,1:end-iN);

%'Estimate B_1,...,B_p'
mBh = adA(mBsh);
mBh =(abs(mBh)>1e-14).*mBh;

%'Estimate A_1,...,A_p'
mB0h = 2*eye(iN) -mC_x(end-iN+1:end,end-iN+1:end);
mB0h =(abs(mB0h)>1e-14).*mB0h;

mAsh = inv(mB0h)*mBsh;
mAh = adA(mAsh);
mAh =(abs(mAh)>1e-14).*mAh;

%'Estimate Sigma_u^{-1}'
mSy = mSx(end-iN+1:end,end-iN+1:end);
mF = [zeros(d,iN) [eye(iN*iP); mAsh]];
mS = zeros(d,d); mS(end-iN+1:end,end-iN+1:end)=mSy;
vSu0 = (eye(d*d)-kron(mF,mF))*vec(mS);
mSu0 = reshape(vSu0,d,d);mSu = mSu0(end-iN+1:end,end-iN+1:end);
mSu =(abs(mSu)>1e-14).*mSu;
%imSu = inv(mSu);
%imSu =(abs(imSu)>1e-14).*imSu;
%mSu = inv(imSu);

end
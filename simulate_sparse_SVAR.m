function [Theta,mAs,mSig_u,mB,mB0] = simulate_sparse_SVAR(N,P,sparsity,sparsity2,density,a1,a2)

% Inputs: - N: dimension of the vector of observations
%        - P: number of desired VAR/SVAR lags
%        - sparsity: desired sparsity degree in Theta_true, i.e., the
%        number corresponding to the proportion of zero coefficients 
%        located in vech(Theta_true), the column vector that stacks the 
%        columns of the lower triangular part of Theta_true, excluding the
%        diagonal coefficients. The proportion of zero entries is simply
%        defined as sparsity/(q*(q-1)/2), q = N(P+1)
%        - sparsity2: controls for the minimum sparsity degree of the zeros
%        in B0, B1, ..., Bp: it refers to \tilde{s} defined step (v) in the
%        simulation procedure described in S.2 of the Online Supplement
%        - density: controls the proportion of non-zero elements in the 
%        q x q matrix, where each entry is the sum of one or more normally 
%        distributed random samples
%        - a1 and a2: coefficients controlling for the generating C^{(1)}
%        in step (ii) of the simulation procedure, i.e. 
%        C^{(1)} = K + \beta I_q, with \beta \in U([a1,a2]) to ensure the
%        positive-definiteness of C^{(0)}, q = N(P+1)

% Outputs: - Theta: q x q matrix precision matrix deduced from the gVAR model
%         with q = N(P+1)
%         - mAs: VAR autoregressive coefficients deduced from Theta
%         - mSig_u: variance-covariance matrix of the error system vector,
%         deduced from B0
%         - mB: sparse SVAR autoregressive matrix coefficients B1,
%         B2,...,Bp
%         - mB0: sparse SVAR matrix B0

% Note: the true output Theta satisfies the following four properties:
%       - positive-definite, such that the VAR process is stationary
%       - sparsity degree in Theta: the proportion of true zero entries in
%       vech(Theta) corresponds to the desired input "sparsity"
%       - sparsity degree in B0, ..., Bp: the poportion of zero entries in
%       B0, ...., Bp should be at least sparsity2 x sparsity/(q*(q-1)/2),
%       with q = N(P+1)
%       - minimum signal condition: the minimum true non-zero coefficient
%       located in vech(Theta), the column vector that stacks the 
%       columns of the lower triangular part of Theta_true, excluding the
%       diagonal coefficients, in absolute value should be larger than 0.05

d = (P+1)*N; cond = true;
while cond
    iE4=0;
    while (iE4==0)
        iE3=0;
        while (iE3==0)
            iE2=0;
            while (iE2==0)
                iE1=0;
                while (iE1==0)
                    mC = full(sprandsym(d,density)+(a1+(a2-a1)*rand(d)).*speye(d));
                    iE1 = (min(eig(mC))>0.1);
                end
                [~,~,mA,~]=Theta2Ba(mC,N);
                mSig_u = inv(mC(end-N+1:end,end-N+1:end));
                mAs = adA(mA);
                mF = [zeros(d,N) [eye(N*P); mAs]];
                iE2 = (max(abs(eig(mF)))<1)&&(min(eig(mSig_u))>0);
            end
            mS = zeros(d,d); mS(end-N+1:end,end-N+1:end)=mSig_u;
            vSig_x = (eye(d*d)-kron(mF,mF))\vec(mS);
            mSig_x = reshape(vSig_x,d,d);
            Theta = inv(mSig_x);
            mR=corrcov(Theta);
            if (P==1)
                Q = {mR(N*P+1:end,N*P+1:end),mR(N*P+1:end,N*(P-1)+1:N*(P-1+1))};
            elseif (P==2)
                Q = {mR(N*P+1:end,N*P+1:end),mR(N*P+1:end,N*(P-1)+1:N*(P-1+1)),mR(N*P+1:end,N*(P-2)+1:N*(P-2+1))};
            elseif (P==3)
                Q = {mR(N*P+1:end,N*P+1:end),mR(N*P+1:end,N*(P-1)+1:N*(P-1+1)),mR(N*P+1:end,N*(P-2)+1:N*(P-2+1)),mR(N*P+1:end,N*(P-3)+1:N*(P-3+1))};
            end
            mQ = cell2mat(Q(toeplitz(1:(P+1))));
            mQ = unvech(vech(mQ));
            vT_q = sort(vech(abs(mQ)));
            iTau_q = vT_q(sparsity);
            mQ = unvech(vech((abs(mQ)>iTau_q).*mQ));
            Theta = (abs(mQ)>0).*Theta;
            vT = sort(vech(abs(Theta)));
            iE3 = (min(eig(Theta))>0.15)&&(vT(sparsity+1)>0);
        end
        [mB,mB0,mA,mSig_u]=Theta2Ba(Theta,N);
        mAs = adA(mA);
        mF = [zeros(d,N) [eye(N*P); mAs]];
        iE4 = (max(abs(eig(mF)))<1);
    end
    K_m = vech(Theta); K_m(K_m==0)=[]; 
    K1 = vec(mB); k1 = length(K1); K1(abs(K1)>0)=[]; 
    K2 = vec(mB0); k2 = length(K2); K2(abs(K2)>0)=[];
    K3 = [K1;K2]; 
    cond = (min(abs(K_m))<0.05)&&((length(K3)/(k1+k2))<sparsity2*(sparsity/(d*(d-1)/2)));
end
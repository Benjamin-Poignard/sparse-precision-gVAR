% Empirical analysis using Realized Kernal (RK)
% mD2: RK for 18 indices: see README and readme.m for a description of the data
% ==> Replication of Table 2 and Figure 1 of the main paper
clear
clc
load('data_rk.mat');
% T: sample size for rolling window
% iF: size of forecasting period
iF = 100;
T = 2000;
% log of RK based on percent (%) return
mYY = log(1e4*mD2(end-T-iF+1:end,2:end));
% information for year, month, and day
vYMD = mD2(end-T-iF+1:end,1);
% mean-sbtracted series
mY0 = mYY -ones(length(mYY),1)*mean(mYY);
[~,N]=size(mY0);

% Set lag length: 1
P=1;

'T,N,P'
[T N P]

% d: dimension of the vector X_t
d = (P+1)*N;
% dim: total number of coefficients in Theta
dim = d^2;
% dim_distinct: total number of distinct coefficients in Theta
dim_distinct = d*(d+1)/2;
dim_distinct2 = d*(d-1)/2;

% To save estimated Theta or A1 matrices
mF_lasso_g = zeros(iF,d*(d+1)/2);
mF_scad_dt = zeros(iF,d*(d+1)/2);
mF_mcp_dt = zeros(iF,d*(d+1)/2);
mF_lasso_dt = zeros(iF,d*(d+1)/2);
mF_alasso_dt = zeros(iF,d*(d+1)/2);
mF_est_lasso = zeros(iF,N^2);
mF_est_ridge = zeros(iF,N^2);

for j=1:iF
    j
    
    % To run the code for real data, the following matrices are necessary
    mY = mY0(j:T+j-1,:);
    
    mX = zeros(T-5,N*(P+1));
    for i=1:(P+1)
        mX(:,(i-1)*N+1:i*N) = mY(i:end-5-1+i,:);
    end
    % mX is the key data input
    Sx = cov(mX); iSx = inv(Sx);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% Gaussian estimator %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The GLASSO algorithm is implemented for Gaussian LASSO precision matrix
    % For computation gain, the grid is set as (0.001:0.01:1) since the
    % GLASSO tends to always select small lambda values

    lambda = (0.001:0.01:0.5)*sqrt(log(dim)/T);
    
    a_scad_g = 3.7; b_mcp_g = 3.5;
    % LASSO
    [theta_est_lasso,lambda_lasso_g_opt]= sparse_precision(mX,'Gaussian',lambda,'lasso',a_scad_g,b_mcp_g);
    % Theta_est_lasso_g: Estimator of Theta based on the Gaussian loss with LASSO penalty
    Theta_est_lasso_g = dvech(theta_est_lasso,N*(P+1));
    % lambda_lasso_g_opt: optimal regularization parameter for LASSO penalty
    % and Gaussian case, identified by cross-validation
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% D-trace squares estimator %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For computation gain, the grid is set as (2.5:0.01:4) since the
    % SCAD/MCP tends to select large lambda values
    
    % calibration of the penalisation parameter for DTrace loss
    lambda = (2.5:0.01:4)*sqrt(log(dim)/T);
    
    % SCAD/MCP coefficients for DTrace loss
    a_scad_dt = 3.7; b_mcp_dt = 3.5;
    
    % SCAD
    [theta_est_scad,lambda_scad_dt_opt]= sparse_precision(mX,'DTrace',lambda,'scad',a_scad_dt,b_mcp_dt);
    % Theta_est_scad_dt: Estimator of Theta based on the DTrace loss with SCAD penalty
    Theta_est_scad_dt = dvech(theta_est_scad,N*(P+1));
    % lambda_scad_dt_opt: optimal regularization parameter for SCAD penalty and
    % DTrace case, identified by cross-validation
    
    % MCP
    [theta_est_mcp,lambda_mcp_dt_opt]= sparse_precision(mX,'DTrace',lambda,'mcp',a_scad_dt,b_mcp_dt);
    % Theta_est_mcp_dt: Estimator of Theta based on the DTrace loss with MCP penalty
    Theta_est_mcp_dt = dvech(theta_est_mcp,N*(P+1));
    % lambda_mcp_dt_opt: optimal regularization parameter for MCP penalty and
    % DTrace case, identified by cross-validation

    lambda = (0.5:0.01:4)*sqrt(log(dim)/T);
    
    % LASSO
    [theta_est_lasso,lambda_lasso_dt_opt]= sparse_precision(mX,'DTrace',lambda,'lasso',a_scad_dt,b_mcp_dt);
    % Theta_est_lasso_dt: Estimator of Theta based on the DTrace loss with LASSO penalty
    Theta_est_lasso_dt = dvech(theta_est_lasso,N*(P+1));
    % lambda_lasso_dt_opt: optimal regularization parameter for LASSO penalty
    % and DTrace case, identified by cross-validation
    
    % Adaptive LASSO
    [theta_est_alasso,lambda_alasso_dt_opt]= sparse_precision(mX,'DTrace',lambda,'alasso',a_scad_dt,b_mcp_dt);
    Theta_est_alasso_dt = dvech(theta_est_alasso,N*(P+1));
    % lambda_alasso_dt_opt: optimal regularization parameter for adaptive LASSO penalty
    
    % VAR estimation with LASSO penalizaton
    lambda = (0.001:0.01:3)*sqrt(log(P*N^2)/T);
    mA_est_lasso= sparse_VAR(mY,P,lambda,'lasso');
    % VAR estimation with Ridge penalization
    lambda = (0.0001:0.01:3)*sqrt(log(P*N^2)/T);
    mA_est_ridge= sparse_VAR(mY,P,lambda,'ridge');
    
    mF_lasso_g(j,:) = vech(Theta_est_lasso_g)';
    mF_scad_dt(j,:) = vech(Theta_est_scad_dt)';
    mF_mcp_dt(j,:) = vech(Theta_est_mcp_dt)';
    mF_lasso_dt(j,:) = vech(Theta_est_lasso_dt)';
    mF_alasso_dt(j,:) = vech(Theta_est_alasso_dt)';
    mF_est_lasso(j,:) = vec(mA_est_lasso)';
    mF_est_ridge(j,:) = vec(mA_est_ridge)';
    
end

save('mR_rk_all_for.mat','mF_lasso_g','mF_scad_dt','mF_mcp_dt','mF_lasso_dt','mF_alasso_dt','mF_est_lasso','mF_est_ridge')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Estimation Sumamary %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gaussian loss based estimators:
% - mF_lasso_g: LASSO penalised

% D-Trace loss based estimators:
% - mF_scad_dt: SCAD penalised
% - mF_mcp_dt: mcp penalised
% - mF_lasso_dt: lasso penalised
% - mF_alasso_dt: adaptive lasso penalised

% OLS estimation of VAR(1) process:
% - mF_est_lasso: LASSO penalised (equation by equation)
% - mF_est_ridge: Ridge regularized (no sparsity is fostered with Ridge)

clear
clc
% names of indices for RK
str ={'S&P 500','FTSE 100','Nikkei 225','DAX','Russel 2000','All Ordinaries','CAC 40','Hang Seng','KOSPI','AEX Index','Swiss Market Index','IBEX 35','S&P CNX Nifty','IPC Mexico','Bovespa Index','S&P/TSX','Euro STOXX 50','FTSE MIB'};
str1 ={'S&P 500(-1)','FTSE 100(-1)','Nikkei 225(-1)','DAX(-1)','Russel 2000(-1)','All Ordinaries(-1)','CAC 40(-1)','Hang Seng(-1)','KOSPI(-1)','AEX Index(-1)','Swiss Market Index(-1)','IBEX 35(-1)','S&P CNX Nifty(-1)','IPC Mexico(-1)','Bovespa Index(-1)','S&P/TSX(-1)','Euro STOXX 50(-1)','FTSE MIB(-1)'};

% mD2: RK for 18 indices 
load('data_rk.mat');
% T: sample size for rolling window
% iF: size of forecasting period
T =2000;
iF = 100;
% log of RK based on percent (%) return
mYY = log(1e4*mD2(end-T-iF+1:end,2:end));
% information for year, month, and day
vYMD = mD2(end-T-iF+1:end,1);
% mean-sbtracted series
mY0 = mYY -ones(length(mYY),1)*mean(mYY);

% Y_t and {Y_{t-1},...,Y_{t-5}}
[iT0,iN] = size(mY0);
mX0 = mY0(6-1:end-1,:);
for i=1:5
    mX0(:,(i-1)*iN+1:i*iN) = mY0(6-i:end-i,:);
end
mY0 = mY0(6:end,:);iT=iT0-5-iF;

% to save MSFE and degree of sparsity
mMSFE = zeros(iF,8);vPa = zeros(8,1);
% P: lag length
P=1;
% load estimates of Theta and A1 matrices
% "mR_rk_all_for.mat" was obtained by "main_rk_all_for.m"
load('mR_rk_all_for.mat'); 
% Save the ratio of times when the coeffcients are non-zero in iF times: 
mB0a = zeros(iN,iN);
mB1a = zeros(iN,iN);
% use estimates of Theta to obtain A_1
% use A_1 to one-step-ahead forecast, Y_{T+1}
for i=1:5
    if (i==1)
        mF = mF_lasso_g;
    elseif (i==2)
        mF = mF_scad_dt;
    elseif (i==3)
        mF = mF_mcp_dt;
    elseif (i==4)
        mF = mF_lasso_dt;
    elseif (i==5)
        mF = mF_alasso_dt;
    end
    vMSE = zeros(iF,1);vQ = zeros(iF,1); vP = zeros(iF,1);
    for k=1:iF
        vTheta = mF(k,:)';mTheta = unvech(vTheta);
        vY = mY0(iT+k,:);
        vY1 = mX0(iT+k,1:P*iN);
        [mBh,mB0h,mAh,mSu]=Theta2Ba(mTheta,iN);
        vYh = vY1*mAh';
        vMSE(k) = sum((vY-vYh).^2);
        [iM,iR]=size(mBh);
        vP(k)=sum(vech_on(abs(mTheta)<1e-14,iM+iR))/(0.5*(iM+iR)*(iM+iR-1));
        if (i==2) % for SCAD D-trace to reproduce Figure 1
            mB0a = mB0a + (abs(mB0h)>0);
            mB1a = mB1a + (abs(mBh)>0);
        end
    end
    mMSFE(:,i)=vMSE;
    vPa(i)=mean(vP);
end
mB0a = (1/iF)*mB0a;
mB1a = (1/iF)*mB1a;

% Plot heatmaps for non-zero coefficients, SCAD-based estimator (i==2 in line 64)
subplot(1,2,1);heatmap(str,str,mB0a);title('(a) Contemporaneous Dependence');
subplot(1,2,2);heatmap(str1,str,mB1a);title('(b) Temporal Dependence')

% save MFSE for OLS, RSS-LASSO, RSS-Ridge
vMSE = zeros(iF,1); vMSE6 = zeros(iF,1); vMSE7 = zeros(iF,1);
for k=1:iF
    mY = mY0(k:iT+k-1,:);
    mY1 = mX0(k:iT+k-1,1:P*iN);
    mAr = (inv(mY1'*mY1)*(mY1'*mY))';
    mE = mY - mY1*mAr';
    imSu = inv((1/iT)*(mE'*mE));iQq=length(vech(imSu));
    size(imSu)
    vY = mY0(iT+k,:);
    vY1 = mX0(iT+k,1:P*iN);
    vYh = vY1*mAr';
    vMSE(k) = sum((vY-vYh).^2);

    mA2 = reshape(mF_est_lasso(k,:),iN,iN);
    vYh = vY1*mA2';
    vMSE6(k) = sum((vY-vYh).^2);
    mA3 = reshape(mF_est_ridge(k,:),iN,iN);
    vYh = vY1*mA3';
    vMSE7(k) = sum((vY-vYh).^2);
end

mMSFE(:,6)=vMSE;mMSFE(:,7)=vMSE6;mMSFE(:,8)=vMSE7;
% the sample means of degree of sparsity and MSFE
format short
% vPa: provides column 3 "Sparsity" of Table 2; how to read this column:
% line 1: GLASSO-based sparsity
% line 2: D-trace SCAD-based sparsity
% line 3: D-trace MCP-based sparsity
% line 4: D-trace LASSO-based sparsity
% line 5: D-trace aLASSO-based sparsity
vPa

vMSFEa = mean(mMSFE)';
% vMSFEa: provides column 4 "MSFE" of Table 2; how to read this column:
% line 1: GLASSO-based MSFE
% line 2: D-trace SCAD-based MSFE
% line 3: D-trace MCP-based MSFE
% line 4: D-trace LASSO-based MSFE
% line 5: D-trace aLASSO-based MSFE
% line 6: VAR OLS-based MSFE
% line 7: LASSO VAR OLS-based MSFE
% line 8: Ridge VAR OLS-based MSFE
vMSFEa

% MCS test: includedR and pvalsR are based on the range statistics provided
% by Hansen, Lunde and Nason
% a bootstrap method is applied, hence set the seed rng(1, 'twister' );
rng(1, 'twister' );
[includedR, pvalsR,excludedR,includedSQ,pvalsSQ,excludedSQ] = mcs(mMSFE,0.05,100000,12);
MCS = [ [excludedR;includedR] , pvalsR ];
% index of each model: 
% 1: GLASSO
% 2: SCAD D-trace
% 3: MCP D-trace
% 4: LASSO
% 5: aLASSO
% 6: OLS VAR
% 7: LASSO OLS VAR
% 8: Ridge OLS VAR

% pvalsR: provides column 5 "MCS" of Table 2: how to read 
% [ [excludedR;includedR] , pvalsR ]:
% - first column: index corresponding to the estimation method ==> refer to
% lines 258-265 above
% - second column: corresponding p-values

clearvars -except vPa vMSFEa MCS
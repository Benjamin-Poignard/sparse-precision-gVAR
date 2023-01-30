% Simulation experiment for N = 10, P = 2 and N = 20, P = 1
% ==> Replication of Table 1 of the main paper
% simulation setting: 90% of the coefficients located below the main
% diagonal of the true precision matrix Theta_true are zero coefficients
% Two sample sizes: T = 500 and T = 2000

% Below are the cases N = 10, p = 2, T = 500 and T = 2000
% N = 10, P = 2, 90% sparse case, sample size T = 500
clear
clc
rng(1, 'twister' );
% N: dimension of the gVAR vector; P: number of lags
N = 10; P = 2;

% d: dimension of the vector X_t; dim: total number of coefficients in the
% precision matrix
d = (P+1)*N; dim = d^2;

% Nsim: number of batches
Nsim = 200;
% T: sample size
T = 500+P;
check = zeros(Nsim,5); check_prop = zeros(Nsim,5);
check_prop2 = zeros(Nsim,5); check2 = zeros(Nsim,10);
loss1 = zeros(Nsim,5); loss2 = zeros(Nsim,5); loss3 = zeros(Nsim,5);

% dim_distinct: number of distinct parameters in the precision matrix (lower triangular part)
dim_distinct = d*(d-1)/2;
% sparsity: number of zero coefficients in dim_distinct, obtained as the
% proportion of true zero coefficients in dim_distinct
sparsity = round(0.9*dim_distinct);
% Nonzero: number of non-zero coefficients in dim_distinct
Nonzero = dim_distinct-sparsity;
% a_1 and a_2 ensure the positive-definiteness of the true precision matrix
a_1 = 4; a_2 = 8;

for oo = 1:Nsim
    
    % simulate a true sparse precision matrix Theta from which we deduce
    % the sparse gVAR parameters according to the supplementary file
    [Theta_true,mAs_true,mSig_u_true] = simulate_sparse_SVAR(N,P,sparsity,0.95,0.2,a_1,a_2);
    
    % generate a sample from the VAR-based observations
    mY0=zeros(T+50,N);
    for t=1:T+50
        vX = vec(mY0(t:t+P-1,:)')';
        mY0(t+P,:) = vX*mAs_true' + random('Normal',0,1,1,N)*sqrtm(mSig_u_true);
    end
    mY = mY0(end-T+1:end,:);
    
    mX = zeros(T-P,N*(P+1));
    for i=1:(P+1)
        mX(:,(i-1)*N+1:i*N) = mY(i:end-P-1+i,:);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% Gaussian estimator %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calibration of the penalisation parameter for estimating a sparse
    % precision matrix
    lambda = (0.001:0.01:3)*sqrt(log(dim)/T);
    
    a_scad_g = 3.7; b_mcp_g = 3.5;
    
    % GLASSO
    [theta_est_lasso,lambda_opt]= sparse_precision(mX,'Gaussian',lambda,'lasso',a_scad_g,b_mcp_g);
    Theta_est_lasso_g = dvech(theta_est_lasso,N*(P+1));
    [N1_lasso_g,N2_lasso_g] = check_IC_2(vech_on(Theta_true,d),vech_on(Theta_est_lasso_g,d));
    NZ_lasso_g = check_NZ(vech_on(Theta_true,d),vech_on(Theta_est_lasso_g,d));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% D-trace squares estimator %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lambda = (0.001:0.01:3)*sqrt(log(dim)/T);
    
    a_scad_dt = 3.7; b_mcp_dt = 3.5;
    
    % SCAD
    [theta_est_scad,~]= sparse_precision(mX,'DTrace',lambda,'scad',a_scad_dt,b_mcp_dt);
    Theta_est_scad_dt = dvech(theta_est_scad,N*(P+1));
    [N1_scad_dt,N2_scad_dt] = check_IC_2(vech_on(Theta_true,d),vech_on(Theta_est_scad_dt,d));
    NZ_scad_dt = check_NZ(vech_on(Theta_true,d),vech_on(Theta_est_scad_dt,d));
    
    % MCP
    [theta_est_mcp,~]= sparse_precision(mX,'DTrace',lambda,'mcp',a_scad_dt,b_mcp_dt);
    Theta_est_mcp_dt = dvech(theta_est_mcp,N*(P+1));
    [N1_mcp_dt,N2_mcp_dt] = check_IC_2(vech_on(Theta_true,d),vech_on(Theta_est_mcp_dt,d));
    NZ_mcp_dt = check_NZ(vech_on(Theta_true,d),vech_on(Theta_est_mcp_dt,d));
    
    % ALASSO
    [theta_est_alasso,~]= sparse_precision(mX,'DTrace',lambda,'alasso',a_scad_dt,b_mcp_dt);
    Theta_est_alasso_dt = dvech(theta_est_alasso,N*(P+1));
    [N1_alasso_dt,N2_alasso_dt] = check_IC_2(vech_on(Theta_true,d),vech_on(Theta_est_alasso_dt,d));
    NZ_alasso_dt = check_NZ(vech_on(Theta_true,d),vech_on(Theta_est_alasso_dt,d));
    
    % LASSO
    [theta_est_lasso,~]= sparse_precision(mX,'DTrace',lambda,'lasso',a_scad_dt,b_mcp_dt);
    Theta_est_lasso_dt = dvech(theta_est_lasso,N*(P+1));
    [N1_lasso_dt,N2_lasso_dt] = check_IC_2(vech_on(Theta_true,d),vech_on(Theta_est_lasso_dt,d));
    NZ_lasso_dt = check_NZ(vech_on(Theta_true,d),vech_on(Theta_est_lasso_dt,d));
    
    % count the number of zero entries correctly identified
    check(oo,:) = [  check_zero(vech_on(Theta_true,d),vech_on(Theta_est_lasso_g,d))...
        check_zero(vech_on(Theta_true,d),vech_on(Theta_est_scad_dt,d)) check_zero(vech_on(Theta_true,d),vech_on(Theta_est_mcp_dt,d))...
        check_zero(vech_on(Theta_true,d),vech_on(Theta_est_alasso_dt,d)) check_zero(vech_on(Theta_true,d),vech_on(Theta_est_lasso_dt,d))];
    
    % count the number of zero entries incorrectly identified and wrongly
    % estimated zero coefficients
    check2(oo,:) = [
        N1_lasso_g N2_lasso_g ...
        N1_scad_dt N2_scad_dt N1_mcp_dt N2_mcp_dt N1_alasso_dt N2_alasso_dt N1_lasso_dt N2_lasso_dt];
    
    % proportion of zero entries correctly estimated
    check_prop(oo,:) = check(oo,:)./sparsity;
    
    % proportion of non-zero entries correctly estimated
    check_prop2(oo,:) = [ NZ_lasso_g  NZ_scad_dt NZ_mcp_dt NZ_alasso_dt NZ_lasso_dt]./Nonzero;
    
    % compute three losses for estimation accuracy: quadratic loss, entropy
    % loss and frobenius norm
    loss1(oo,:) = [quadratic_loss(Theta_true,Theta_est_lasso_g)...
        quadratic_loss(Theta_true,Theta_est_scad_dt) quadratic_loss(Theta_true,Theta_est_mcp_dt) ...
        quadratic_loss(Theta_true,Theta_est_alasso_dt) quadratic_loss(Theta_true,Theta_est_lasso_dt) ];
    
    
    loss2(oo,:) = [ entropy_loss(Theta_true,Theta_est_lasso_g)...
        entropy_loss(Theta_true,Theta_est_scad_dt) entropy_loss(Theta_true,Theta_est_mcp_dt) ...
        entropy_loss(Theta_true,Theta_est_alasso_dt) entropy_loss(Theta_true,Theta_est_lasso_dt) ];
    
    loss3(oo,:) = [ norm(Theta_true-Theta_est_lasso_g,'fro')...
        norm(Theta_true-Theta_est_scad_dt,'fro') norm(Theta_true-Theta_est_mcp_dt,'fro') ...
        norm(Theta_true-Theta_est_alasso_dt,'fro') norm(Theta_true-Theta_est_lasso_dt,'fro') ];
       
end

% Report the result of Table 1, case T = 500, N = 10, p = 2
% Results_1 has lines and columns
%              GLASSO    D-T SCAD    D-T MCP    D-T aLASSO    D-T LASSO
% mean(loss1)
% std(loss1)
% mean(loss2)
% std(loss2)
% mean(loss3)
% std(loss3)
% mean(C1)
% std(C1)
% mean(C2)
% std(C2)

% Remark: the results in Table 1 of the main paper, for each column, are 
% displayed with the mean and standard deviation in parenthesis for each
% line
Results_1 = round([
    mean(loss1);std(loss1);mean(loss2);std(loss2);mean(loss3);std(loss3);...
    mean(100*check_prop);std(100*check_prop);...
    mean(100*check_prop2);std(100*check_prop2);...
    ],3);

% N = 10, P = 2, 90% sparse case, sample size T = 2000
clearvars -except Results_1
rng(2, 'twister' );
N = 10; P = 2;

d = (P+1)*N; dim = d^2;

Nsim = 200;
T = 2000+P;
check = zeros(Nsim,5); check_prop = zeros(Nsim,5);
check_prop2 = zeros(Nsim,5); check2 = zeros(Nsim,10);
loss1 = zeros(Nsim,5); loss2 = zeros(Nsim,5); loss3 = zeros(Nsim,5);
dim_distinct = d*(d-1)/2; sparsity = round(0.9*dim_distinct);
Nonzero = dim_distinct-sparsity;
a_1 = 4; a_2 = 8;

for oo = 1:Nsim
    
    [Theta_true,mAs_true,mSig_u_true] = simulate_sparse_SVAR(N,P,sparsity,0.95,0.2,a_1,a_2);
    mY0=zeros(T+50,N);
    for t=1:T+50
        vX = vec(mY0(t:t+P-1,:)')';
        mY0(t+P,:) = vX*mAs_true' + random('Normal',0,1,1,N)*sqrtm(mSig_u_true);
    end
    mY = mY0(end-T+1:end,:);
    
    mX = zeros(T-P,N*(P+1));
    for i=1:(P+1)
        mX(:,(i-1)*N+1:i*N) = mY(i:end-P-1+i,:);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% Gaussian estimator %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calibration of the penalisation parameter for estimating a sparse
    % precision matrix
    lambda = (0.001:0.01:3)*sqrt(log(dim)/T);
    
    a_scad_g = 3.7; b_mcp_g = 3.5;
    
    % GLASSO
    [theta_est_lasso,lambda_opt]= sparse_precision(mX,'Gaussian',lambda,'lasso',a_scad_g,b_mcp_g);
    Theta_est_lasso_g = dvech(theta_est_lasso,N*(P+1));
    [N1_lasso_g,N2_lasso_g] = check_IC_2(vech_on(Theta_true,d),vech_on(Theta_est_lasso_g,d));
    NZ_lasso_g = check_NZ(vech_on(Theta_true,d),vech_on(Theta_est_lasso_g,d));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% D-trace squares estimator %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lambda = (0.001:0.01:3)*sqrt(log(dim)/T);
    
    a_scad_dt = 3.7; b_mcp_dt = 3.5;
    
    % SCAD
    [theta_est_scad,~]= sparse_precision(mX,'DTrace',lambda,'scad',a_scad_dt,b_mcp_dt);
    Theta_est_scad_dt = dvech(theta_est_scad,N*(P+1));
    [N1_scad_dt,N2_scad_dt] = check_IC_2(vech_on(Theta_true,d),vech_on(Theta_est_scad_dt,d));
    NZ_scad_dt = check_NZ(vech_on(Theta_true,d),vech_on(Theta_est_scad_dt,d));
    
    % MCP
    [theta_est_mcp,~]= sparse_precision(mX,'DTrace',lambda,'mcp',a_scad_dt,b_mcp_dt);
    Theta_est_mcp_dt = dvech(theta_est_mcp,N*(P+1));
    [N1_mcp_dt,N2_mcp_dt] = check_IC_2(vech_on(Theta_true,d),vech_on(Theta_est_mcp_dt,d));
    NZ_mcp_dt = check_NZ(vech_on(Theta_true,d),vech_on(Theta_est_mcp_dt,d));
    
    % ALASSO
    [theta_est_alasso,~]= sparse_precision(mX,'DTrace',lambda,'alasso',a_scad_dt,b_mcp_dt);
    Theta_est_alasso_dt = dvech(theta_est_alasso,N*(P+1));
    [N1_alasso_dt,N2_alasso_dt] = check_IC_2(vech_on(Theta_true,d),vech_on(Theta_est_alasso_dt,d));
    NZ_alasso_dt = check_NZ(vech_on(Theta_true,d),vech_on(Theta_est_alasso_dt,d));
    
    % LASSO
    [theta_est_lasso,~]= sparse_precision(mX,'DTrace',lambda,'lasso',a_scad_dt,b_mcp_dt);
    Theta_est_lasso_dt = dvech(theta_est_lasso,N*(P+1));
    [N1_lasso_dt,N2_lasso_dt] = check_IC_2(vech_on(Theta_true,d),vech_on(Theta_est_lasso_dt,d));
    NZ_lasso_dt = check_NZ(vech_on(Theta_true,d),vech_on(Theta_est_lasso_dt,d));
    
    check(oo,:) = [  check_zero(vech_on(Theta_true,d),vech_on(Theta_est_lasso_g,d))...
        check_zero(vech_on(Theta_true,d),vech_on(Theta_est_scad_dt,d)) check_zero(vech_on(Theta_true,d),vech_on(Theta_est_mcp_dt,d))...
        check_zero(vech_on(Theta_true,d),vech_on(Theta_est_alasso_dt,d)) check_zero(vech_on(Theta_true,d),vech_on(Theta_est_lasso_dt,d))];
    
    check2(oo,:) = [
        N1_lasso_g N2_lasso_g ...
        N1_scad_dt N2_scad_dt N1_mcp_dt N2_mcp_dt N1_alasso_dt N2_alasso_dt N1_lasso_dt N2_lasso_dt];
    
    check_prop(oo,:) = check(oo,:)./sparsity;
    
    check_prop2(oo,:) = [ NZ_lasso_g  NZ_scad_dt NZ_mcp_dt NZ_alasso_dt NZ_lasso_dt]./Nonzero;
    
    loss1(oo,:) = [quadratic_loss(Theta_true,Theta_est_lasso_g)...
        quadratic_loss(Theta_true,Theta_est_scad_dt) quadratic_loss(Theta_true,Theta_est_mcp_dt) ...
        quadratic_loss(Theta_true,Theta_est_alasso_dt) quadratic_loss(Theta_true,Theta_est_lasso_dt) ];
    
    loss2(oo,:) = [ entropy_loss(Theta_true,Theta_est_lasso_g)...
        entropy_loss(Theta_true,Theta_est_scad_dt) entropy_loss(Theta_true,Theta_est_mcp_dt) ...
        entropy_loss(Theta_true,Theta_est_alasso_dt) entropy_loss(Theta_true,Theta_est_lasso_dt) ];
    
    loss3(oo,:) = [ norm(Theta_true-Theta_est_lasso_g,'fro')...
        norm(Theta_true-Theta_est_scad_dt,'fro') norm(Theta_true-Theta_est_mcp_dt,'fro') ...
        norm(Theta_true-Theta_est_alasso_dt,'fro') norm(Theta_true-Theta_est_lasso_dt,'fro') ];
    
end

% Report the result of Table 1, case T = 2000, N = 10, p = 2
% Results_2 has lines and columns
%              GLASSO    D-T SCAD    D-T MCP    D-T aLASSO    D-T LASSO
% mean(loss1)
% std(loss1)
% mean(loss2)
% std(loss2)
% mean(loss3)
% std(loss3)
% mean(C1)
% std(C1)
% mean(C2)
% std(C2)
Results_2 = round([
    mean(loss1);std(loss1);mean(loss2);std(loss2);mean(loss3);std(loss3);...
    mean(100*check_prop);std(100*check_prop);...
    mean(100*check_prop2);std(100*check_prop2);...
    ],3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below are the cases N = 20, p = 1, T = 500 and T = 2000
% N = 20, P = 1, 90% sparse case, sample size T = 500
clearvars -except Results_1 Results_2
rng(3, 'twister' );
% N: dimension of the gVAR vector; P: number of lags
N = 20; P = 1;

% d: dimension of the vector X_t; dim: total number of coefficients in the
% precision matrix
d = (P+1)*N; dim = d^2;

% Nsim: number of batches
Nsim = 200;
% T: sample size
T = 500+P;
check = zeros(Nsim,5); check_prop = zeros(Nsim,5);
check_prop2 = zeros(Nsim,5); check2 = zeros(Nsim,10);
loss1 = zeros(Nsim,5); loss2 = zeros(Nsim,5); loss3 = zeros(Nsim,5);

% dim_distinct: number of distinct parameters in the precision matrix (lower triangular part)
dim_distinct = d*(d-1)/2;
% sparsity: number of zero coefficients in dim_distinct, obtained as the
% proportion of true zero coefficients in dim_distinct
sparsity = round(0.9*dim_distinct);
% Nonzero: number of non-zero coefficients in dim_distinct
Nonzero = dim_distinct-sparsity;
% a_1 and a_2 ensure the positive-definiteness of the true precision matrix
a_1 = 5; a_2 = 9;

for oo = 1:Nsim
    
    % simulate a true sparse precision matrix Theta from which we deduce
    % the sparse gVAR parameters according to the supplementary file
    [Theta_true,mAs_true,mSig_u_true] = simulate_sparse_SVAR(N,P,sparsity,0.95,0.2,a_1,a_2);
    
    % generate a sample from the VAR-based observations
    mY0=zeros(T+50,N);
    for t=1:T+50
        vX = vec(mY0(t:t+P-1,:)')';
        mY0(t+P,:) = vX*mAs_true' + random('Normal',0,1,1,N)*sqrtm(mSig_u_true);
    end
    mY = mY0(end-T+1:end,:);
    
    mX = zeros(T-P,N*(P+1));
    for i=1:(P+1)
        mX(:,(i-1)*N+1:i*N) = mY(i:end-P-1+i,:);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% Gaussian estimator %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calibration of the penalisation parameter for estimating a sparse
    % precision matrix
    lambda = (0.001:0.01:3)*sqrt(log(dim)/T);
    
    a_scad_g = 3.7; b_mcp_g = 3.5;
    
    % GLASSO
    [theta_est_lasso,lambda_opt]= sparse_precision(mX,'Gaussian',lambda,'lasso',a_scad_g,b_mcp_g);
    Theta_est_lasso_g = dvech(theta_est_lasso,N*(P+1));
    [N1_lasso_g,N2_lasso_g] = check_IC_2(vech_on(Theta_true,d),vech_on(Theta_est_lasso_g,d));
    NZ_lasso_g = check_NZ(vech_on(Theta_true,d),vech_on(Theta_est_lasso_g,d));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% D-trace squares estimator %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lambda = (0.001:0.01:3)*sqrt(log(dim)/T);
    
    a_scad_dt = 3.7; b_mcp_dt = 3.5;
    
    % SCAD
    [theta_est_scad,~]= sparse_precision(mX,'DTrace',lambda,'scad',a_scad_dt,b_mcp_dt);
    Theta_est_scad_dt = dvech(theta_est_scad,N*(P+1));
    [N1_scad_dt,N2_scad_dt] = check_IC_2(vech_on(Theta_true,d),vech_on(Theta_est_scad_dt,d));
    NZ_scad_dt = check_NZ(vech_on(Theta_true,d),vech_on(Theta_est_scad_dt,d));
    
    % MCP
    [theta_est_mcp,~]= sparse_precision(mX,'DTrace',lambda,'mcp',a_scad_dt,b_mcp_dt);
    Theta_est_mcp_dt = dvech(theta_est_mcp,N*(P+1));
    [N1_mcp_dt,N2_mcp_dt] = check_IC_2(vech_on(Theta_true,d),vech_on(Theta_est_mcp_dt,d));
    NZ_mcp_dt = check_NZ(vech_on(Theta_true,d),vech_on(Theta_est_mcp_dt,d));
    
    % ALASSO
    [theta_est_alasso,~]= sparse_precision(mX,'DTrace',lambda,'alasso',a_scad_dt,b_mcp_dt);
    Theta_est_alasso_dt = dvech(theta_est_alasso,N*(P+1));
    [N1_alasso_dt,N2_alasso_dt] = check_IC_2(vech_on(Theta_true,d),vech_on(Theta_est_alasso_dt,d));
    NZ_alasso_dt = check_NZ(vech_on(Theta_true,d),vech_on(Theta_est_alasso_dt,d));
    
    % LASSO
    [theta_est_lasso,~]= sparse_precision(mX,'DTrace',lambda,'lasso',a_scad_dt,b_mcp_dt);
    Theta_est_lasso_dt = dvech(theta_est_lasso,N*(P+1));
    [N1_lasso_dt,N2_lasso_dt] = check_IC_2(vech_on(Theta_true,d),vech_on(Theta_est_lasso_dt,d));
    NZ_lasso_dt = check_NZ(vech_on(Theta_true,d),vech_on(Theta_est_lasso_dt,d));
    
    % count the number of zero entries correctly identified
    check(oo,:) = [  check_zero(vech_on(Theta_true,d),vech_on(Theta_est_lasso_g,d))...
        check_zero(vech_on(Theta_true,d),vech_on(Theta_est_scad_dt,d)) check_zero(vech_on(Theta_true,d),vech_on(Theta_est_mcp_dt,d))...
        check_zero(vech_on(Theta_true,d),vech_on(Theta_est_alasso_dt,d)) check_zero(vech_on(Theta_true,d),vech_on(Theta_est_lasso_dt,d))];
    
    % count the number of zero entries incorrectly identified and wrongly
    % estimated zero coefficients
    check2(oo,:) = [
        N1_lasso_g N2_lasso_g ...
        N1_scad_dt N2_scad_dt N1_mcp_dt N2_mcp_dt N1_alasso_dt N2_alasso_dt N1_lasso_dt N2_lasso_dt];
    
    % proportion of zero entries correctly estimated
    check_prop(oo,:) = check(oo,:)./sparsity;
    
    % proportion of non-zero entries correctly estimated
    check_prop2(oo,:) = [ NZ_lasso_g  NZ_scad_dt NZ_mcp_dt NZ_alasso_dt NZ_lasso_dt]./Nonzero;
    
    
    % compute three losses for estimation accuracy: quadratic loss, entropy
    % loss and frobenius norm
    loss1(oo,:) = [quadratic_loss(Theta_true,Theta_est_lasso_g)...
        quadratic_loss(Theta_true,Theta_est_scad_dt) quadratic_loss(Theta_true,Theta_est_mcp_dt) ...
        quadratic_loss(Theta_true,Theta_est_alasso_dt) quadratic_loss(Theta_true,Theta_est_lasso_dt) ];
    
    loss2(oo,:) = [ entropy_loss(Theta_true,Theta_est_lasso_g)...
        entropy_loss(Theta_true,Theta_est_scad_dt) entropy_loss(Theta_true,Theta_est_mcp_dt) ...
        entropy_loss(Theta_true,Theta_est_alasso_dt) entropy_loss(Theta_true,Theta_est_lasso_dt) ];
    
    loss3(oo,:) = [ norm(Theta_true-Theta_est_lasso_g,'fro')...
        norm(Theta_true-Theta_est_scad_dt,'fro') norm(Theta_true-Theta_est_mcp_dt,'fro') ...
        norm(Theta_true-Theta_est_alasso_dt,'fro') norm(Theta_true-Theta_est_lasso_dt,'fro') ];
    
end

% Report the result of Table 1, case T = 500, N = 20, p = 1
% Results_3 has lines and columns 
%              GLASSO    D-T SCAD    D-T MCP    D-T aLASSO    D-T LASSO
% mean(loss1)
% std(loss1)
% mean(loss2)
% std(loss2)
% mean(loss3)
% std(loss3)
% mean(C1)
% std(C1)
% mean(C2)
% std(C2)
Results_3 = round([
    mean(loss1);std(loss1);mean(loss2);std(loss2);mean(loss3);std(loss3);...
    mean(100*check_prop);std(100*check_prop);...
    mean(100*check_prop2);std(100*check_prop2);...
    ],3);

% N = 20, P = 1, 90% sparse case, sample size T = 2000
clearvars -except Results_1 Results_2 Results_3
clc
rng(4, 'twister' );
N = 20; P = 1;

d = (P+1)*N; dim = d^2;

Nsim = 200;
T = 2000+P;
check = zeros(Nsim,5); check_prop = zeros(Nsim,5);
check_prop2 = zeros(Nsim,5); check2 = zeros(Nsim,10);
loss1 = zeros(Nsim,5); loss2 = zeros(Nsim,5); loss3 = zeros(Nsim,5);
dim_distinct = d*(d-1)/2; sparsity = round(0.9*dim_distinct);
Nonzero = dim_distinct-sparsity;
a_1 = 5; a_2 = 9;

for oo = 1:Nsim
    
    [Theta_true,mAs_true,mSig_u_true] = simulate_sparse_SVAR(N,P,sparsity,0.95,0.2,a_1,a_2);
    mY0=zeros(T+50,N);
    for t=1:T+50
        vX = vec(mY0(t:t+P-1,:)')';
        mY0(t+P,:) = vX*mAs_true' + random('Normal',0,1,1,N)*sqrtm(mSig_u_true);
    end
    mY = mY0(end-T+1:end,:);
    
    mX = zeros(T-P,N*(P+1));
    for i=1:(P+1)
        mX(:,(i-1)*N+1:i*N) = mY(i:end-P-1+i,:);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% Gaussian estimator %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calibration of the penalisation parameter for estimating a sparse
    % precision matrix
    lambda = (0.001:0.01:3)*sqrt(log(dim)/T);
    
    a_scad_g = 3.7; b_mcp_g = 3.5;
    
    % GLASSO
    [theta_est_lasso,lambda_opt]= sparse_precision(mX,'Gaussian',lambda,'lasso',a_scad_g,b_mcp_g);
    Theta_est_lasso_g = dvech(theta_est_lasso,N*(P+1));
    [N1_lasso_g,N2_lasso_g] = check_IC_2(vech_on(Theta_true,d),vech_on(Theta_est_lasso_g,d));
    NZ_lasso_g = check_NZ(vech_on(Theta_true,d),vech_on(Theta_est_lasso_g,d));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% D-trace squares estimator %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lambda = (0.001:0.01:3)*sqrt(log(dim)/T);
    
    a_scad_dt = 3.7; b_mcp_dt = 3.5;
    
    % SCAD
    [theta_est_scad,~]= sparse_precision(mX,'DTrace',lambda,'scad',a_scad_dt,b_mcp_dt);
    Theta_est_scad_dt = dvech(theta_est_scad,N*(P+1));
    [N1_scad_dt,N2_scad_dt] = check_IC_2(vech_on(Theta_true,d),vech_on(Theta_est_scad_dt,d));
    NZ_scad_dt = check_NZ(vech_on(Theta_true,d),vech_on(Theta_est_scad_dt,d));
    
    % MCP
    [theta_est_mcp,~]= sparse_precision(mX,'DTrace',lambda,'mcp',a_scad_dt,b_mcp_dt);
    Theta_est_mcp_dt = dvech(theta_est_mcp,N*(P+1));
    [N1_mcp_dt,N2_mcp_dt] = check_IC_2(vech_on(Theta_true,d),vech_on(Theta_est_mcp_dt,d));
    NZ_mcp_dt = check_NZ(vech_on(Theta_true,d),vech_on(Theta_est_mcp_dt,d));
    
    % ALASSO
    [theta_est_alasso,~]= sparse_precision(mX,'DTrace',lambda,'alasso',a_scad_dt,b_mcp_dt);
    Theta_est_alasso_dt = dvech(theta_est_alasso,N*(P+1));
    [N1_alasso_dt,N2_alasso_dt] = check_IC_2(vech_on(Theta_true,d),vech_on(Theta_est_alasso_dt,d));
    NZ_alasso_dt = check_NZ(vech_on(Theta_true,d),vech_on(Theta_est_alasso_dt,d));
    
    % LASSO
    [theta_est_lasso,~]= sparse_precision(mX,'DTrace',lambda,'lasso',a_scad_dt,b_mcp_dt);
    Theta_est_lasso_dt = dvech(theta_est_lasso,N*(P+1));
    [N1_lasso_dt,N2_lasso_dt] = check_IC_2(vech_on(Theta_true,d),vech_on(Theta_est_lasso_dt,d));
    NZ_lasso_dt = check_NZ(vech_on(Theta_true,d),vech_on(Theta_est_lasso_dt,d));
    
    check(oo,:) = [  check_zero(vech_on(Theta_true,d),vech_on(Theta_est_lasso_g,d))...
        check_zero(vech_on(Theta_true,d),vech_on(Theta_est_scad_dt,d)) check_zero(vech_on(Theta_true,d),vech_on(Theta_est_mcp_dt,d))...
        check_zero(vech_on(Theta_true,d),vech_on(Theta_est_alasso_dt,d)) check_zero(vech_on(Theta_true,d),vech_on(Theta_est_lasso_dt,d))];
    
    check2(oo,:) = [
        N1_lasso_g N2_lasso_g ...
        N1_scad_dt N2_scad_dt N1_mcp_dt N2_mcp_dt N1_alasso_dt N2_alasso_dt N1_lasso_dt N2_lasso_dt];
    
    check_prop(oo,:) = check(oo,:)./sparsity;
    
    check_prop2(oo,:) = [ NZ_lasso_g  NZ_scad_dt NZ_mcp_dt NZ_alasso_dt NZ_lasso_dt]./Nonzero;
    
    loss1(oo,:) = [quadratic_loss(Theta_true,Theta_est_lasso_g)...
        quadratic_loss(Theta_true,Theta_est_scad_dt) quadratic_loss(Theta_true,Theta_est_mcp_dt) ...
        quadratic_loss(Theta_true,Theta_est_alasso_dt) quadratic_loss(Theta_true,Theta_est_lasso_dt) ];
    
    loss2(oo,:) = [ entropy_loss(Theta_true,Theta_est_lasso_g)...
        entropy_loss(Theta_true,Theta_est_scad_dt) entropy_loss(Theta_true,Theta_est_mcp_dt) ...
        entropy_loss(Theta_true,Theta_est_alasso_dt) entropy_loss(Theta_true,Theta_est_lasso_dt) ];
    
    loss3(oo,:) = [ norm(Theta_true-Theta_est_lasso_g,'fro')...
        norm(Theta_true-Theta_est_scad_dt,'fro') norm(Theta_true-Theta_est_mcp_dt,'fro') ...
        norm(Theta_true-Theta_est_alasso_dt,'fro') norm(Theta_true-Theta_est_lasso_dt,'fro') ];
    
end

% Report the result of Table 1, case T = 2000, N = 20, p = 1
% Results_4 has lines and columns
%              GLASSO    D-T SCAD    D-T MCP    D-T aLASSO    D-T LASSO
% mean(loss1)
% std(loss1)
% mean(loss2)
% std(loss2)
% mean(loss3)
% std(loss3)
% mean(C1)
% std(C1)
% mean(C2)
% std(C2)
Results_4 = round([
    mean(loss1);std(loss1);mean(loss2);std(loss2);mean(loss3);std(loss3);...
    mean(100*check_prop);std(100*check_prop);...
    mean(100*check_prop2);std(100*check_prop2);...
    ],3);

clearvars -except Results_1 Results_2 Results_3 Results_4
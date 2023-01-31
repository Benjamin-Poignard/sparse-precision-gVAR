
# Estimation Vector Autoregression via Sparse Precision Matrix

Matlab implementation of high-dimensional Vector Autoregression based on sparse precision matrix replicating the results of the paper:

Estimation of High-Dimensional Vector Autoregression via Sparse Precision Matrix, by Benjamin Poignard and Manabu Asai.

Link: https://doi.org/10.1093/ectj/utad003

# Overview

The code in this replication package reproduces:

- The results provided in Table 1 based on the simulated experiment detailed in Section 4: the replicator should execute program *simulations.m*.  
    
- The results provided in Table 2 and Figure 1 (provided in the replication package) based on the real data analysis detailed in Section 5: the replication should execute program *real_data_analysis.m*.

The replicator can find a brief overview of the programs to run in *readme.m*.


# Data Availability and Provenance Statements

The data used to support the findings of this study are publicly available. The data were collected by the authors from the database "Oxford-Man Institute's realized library", Version 0.3, produced by Heber et al. (2009), from the link (free access, given "as is"): https://realized.oxford-man.ox.ac.uk (accessed: 21 August 2022). The authors certify that they have legitimate access to and permission to use the data.
 
The raw data file to be downloaded from the Oxford Man Institute of Quantitative Finance is *data_rk.csv* (provided in the replication package). The Oxford Man Institute's realized library contains daily non-parametric metrics of past volatility financial indices. 

The realized volatility measure for the index STI is available until Sept. 18, 2015, while the remaining series contain data up to Dec. 4, 2017. To obtain *data_rk.mat* provided in the package, the replicator should run file *create_data_rk.m*: *data_rk.mat* contains two sets of data:

- mD1: All data starting from Jan. 3, 2000 to Sep. 18, 2015.
- mD2: All data, excluding STI, starting from Jan. 3, 2000 to Dec. 4, 2017.

The missing data were fulfilled using the value of the previous day. The real data analysis is based on mD2 for the period from Nov. 23, 2009 to Dec. 4, 2017. The list of indices to be found in *data_rk.csv*, and composing mD2 excluding STI, is as follows: SPX for S&P 500 (USA); FTSE for FTSE 100 (UK); N225 for Nikkei 225 (Japan); GDAXI for DAX (Germany); RUT for Russel 2000 (USA); AORD for All Ordinaries (Australia); FCHI for CAC 40 (France); HSI for Hang Seng (Hong Kong); KS11 for Korea Composite Stock Price Index - KOSPI (South Korea); AEX for AEX Index (The Netherlands); SSMI for Swiss Market Index (Switzerland); IBEX for IBEX 35 (Spain); NSEI for S&P CNX Nifty 50 (India); MXX for IPC Mexico (Mexico); BVSP for Bovespa Index (Brazil); GSPTSE for S&P/TSX (Canada); STOXX50E for Euro STOXX 50 (50 largest eurozone caps); STI for FT Straits Times Index (Singapore); FTMIB for FTSE MIB (Italy).

# Computational Requirements

**Software requirements**

The Matlab code was run on a Mac-OS Apple M1 Ultra with 20 cores and 128 GB Memory. The version of the Matlab software on which the code was run is a follows: 9.12.0.1975300 (R2022a) Update 3. 

The following toolboxes should be installed:
- Statistics and Machine Learning Toolbox, Version 12.3.
- Parallel Computing Toolbox, Version 7.6.
Parallel Computing Toolbox is highly recommended to run the code to speed up the cross-validation procedure employed to select the optimal tuning parameter. All the run-time requirements displayed below are reported when the code is run with the Parallel Computing Toolbox.


**Remark 1.** The Version 12.3 (or Version 12.0) of the Statistics and Machine Learning Toolbox allows to specify no-intercept in the function "lasso.m": this is important to be able to run *real_data_analysis.m* to get the LASSO-OLS VAR estimates. Indeed,  Version 11.6 does not able for the user to specify "no intercept". But the user may directly implement the LASSO solution in the function *VAR_penalized.m* to avoid the use of the Statistics and Machine Learning Toolbox.

**Remark 2.** No Toolbox is required to run the functions *dtrace_admm.m* (D-trace loss method) and *graphicalLasso.m* (GLASSO): see **C. Description of the main functions** for a description. No Toolbox is required to run *simulations.m*, but the Parallel Computing Toolbox is strongly recommended.

**Controlled randomness**

For the sake of replication, seeds have been specified in *simulations.m* and *real_data_analysis.m*. To be specific:
- In *simulations.m*, the seed is set at:
  - line 11 for the case $N=10,p=2,T = 500$.
  - line 161 for the case $N=10,p=2,T = 2000$.
  - line 287 for the case $N=20,p=1,T = 500$.
  - line 434 for the case $N=20,p=1,T = 2000$.

- In *real_data_analysis.m*, the seed is set at line $271$ since the MCS test is based on a bootstrap procedure.


**Memory and Run-time Requirements**

The approximate time needed to reproduce the simulated experiment and the real data experiment on a Mac-OS M1 Ultra desktop machine is as follows:

- To run file *simulations.m* to obtain all the results reported in Table 1, the approximate computation time is: $9$ hours $45$ mins. For each batch, the approximate computation time for each penalization case is as follows:
    
    - when $N=10,p=2,T=500$:

        GLASSO: $16$ to $20$ seconds; SCAD/MCP D-trace: $4$ to $7$ seconds; adaptive LASSO D-trace: $4$ to $6$ seconds; LASSO D-trace: $2$ to $4$ seconds. The approximate run-times for $T=2000$ are similar.
   
   - when $N=20,p=1,T=500$:

        GLASSO: $22$ to $25$ seconds; SCAD/MCP D-trace: $13$ to $16$ seconds; adaptive LASSO D-trace: $8$ to $10$ seconds; LASSO D-trace: $6$ to $8$ seconds. The approximate run-times for $T=2000$ are similar.


- To run file *real_data_analysis.m* to obtain all the results reported in Table 2 and replicate Figure 1, the approximate computation time is: $12$ hours and $15$ minutes. 
The approximate computation time for a one-step ahead estimation and for each penalization case, is as follows:
GLASSO: $5$ to $6$ minutes; SCAD/MCP D-trace: $25$ to $30$ seconds; adaptive LASSO D-trace: $25$ to $30$ seconds; LASSO D-trace: $20$ to $24$ seconds.

**Remark 3.** The computation time highly depends on the grid size selected for cross-validation to choose the optimal $\lambda_T$ (called tuning or regularization parameter), denoted by $\lambda^\ast_T$: in all the simulated experiments, the optimal $\lambda^\ast_T$ is searched in the grid $\{c\sqrt{\log(q)/T}, c=0.001:0.01:3\}$, so that there are $300$ different $\lambda_T$ candidates.
To save run-time computation, the user may set a smaller grid, such as $\{c\sqrt{\log(q)/T}, c=0.001:0.1:3\}$; the upper bound value of the grid can be smaller for the LASSO-based method, since the LASSO tends to select small $\lambda_T$ values. A larger lower bound value of the grid can be set for SCAD/MCP-based methods, as they tend to select larger $\lambda_T$ values: this is why different grid values have been specified for GLASSO, SCAD/MCP D-trace and adaptive LASSO/LASSO D-trace in *real_data_analysis.m*.

**Remark 4.** In both *dtrace_admm.m* and *graphicalLasso.m*, the maximum number of iterations is set as $\text{maxIt} = 100000$, by default. The user may set smaller values directly in both codes. 

For the sake of comparison in terms of run-time computation, the approximate computation time reported below holds when the code was run on a PC Precision 7920 Tower(Xeon Silver4210) with CPU Xeon Silver 4214 (2.2 GHz, 12 Core), with Matlab Version 9.7.0.1296695 (R2019b) Update 4, with toolboxes:

- Statistics and Machine Learning Toolbox, 11.6.
- Parallel Computing Toolbox, Version 7.1.

In light of **Remark 1**, the function *sparse_VAR.m* in *real_data_analysis.m*, line 114, can not be executed for the LASSO case with Version 11.6 of the Statistics and Machine Learning Toolbox. Thus, we report here the computation-time when running *simulations.m* and *real_data_analysis.m*, excluding "mA_est_lasso= sparse_VAR(mY,P,lambda,'lasso')" (line 114) in the latter file.

- To run file *simulations.m* to obtain all the results reported in Table 1, the approximate computation time is: $31$ hours and $10$ minutes. For each batch, the approximate computation time for each penalization case, is as follows:
    
    - when $N=10,p=2,T=500$:
    
      GLASSO: $35$ to $40$ seconds; SCAD/MCP D-trace: $22$ to $26$ seconds; adaptive LASSO D-trace: $12$ to $16$ seconds; LASSO D-trace: $9$ to $12$ seconds. The approximate run-times for $T=2000$ are similar.

    - when $N=20,p=1,T=500$:

        GLASSO: $45$ to $50$ seconds; SCAD/MCP D-trace: $40$ to $44$ seconds; adaptive LASSO D-trace: $23$ to $25$ seconds; LASSO D-trace: $18$ to $20$ seconds. The approximate run-times for $T=2000$ are similar.

- To run file *real_data_analysis.m* to obtain all the results reported in Table 2 and replicate Figure 1, the approximate computation time is: $29$ hours and $25$ minutes. 

The approximate computation time for a one-step ahead estimation and for each penalization case, is as follows:
GLASSO: $9$ minutes to $9$ minutes and $30$ seconds; SCAD/MCP D-trace: $1$ minute and $30$ seconds to $2$ minutes; adaptive LASSO D-trace: $1$ minute and $50$ seconds to $2$ minutes and $30$ seconds; LASSO D-trace: $1$ minute and $10$ to $20$ seconds.


# Description of the code for replication

The replicator should execute the following Matlab files to replicate Table 1, Table 2, Figure 1 of the paper.

**A. Replication of the results reported in Table 1 of the main paper**

Program *simulations.m* will replicate the results reported in Table 1. Specifically, the following matrix variables will be produced:
  - "Results_1" (line $153$) replicates the results of Table 1 for $N=10,p=2,T=500$.
  - "Results_2" (line $275$) replicates the results of Table 1 for $N=10,p=2,T=2000$.
  - "Results_3" (line $425$) replicates the results of Table 1 for $N=20,p=1,T=500$.
  - "Results_4" (line $548$) replicates the results of Table 1 for $N=20,p=1,T=2000$.

**B. Replication of the results reported in Table 2 and Figure 1 of the main paper**

- Program *real_data_analysis.m* will replicate the results reported in Table 2 and produce Figure 1. Specifically, the program performs the rolling window estimation of the gVAR with sparse precision matrix (GLASSO and D-trace with SCAD, MCP, LASSO and adaptive LASSO) and of the VAR with LASSO OLS and Ridge OLS.
- Once the program has been executed, the replicator will find:
  - the vector variable "vPa" (line 214): provides column 3 "Sparsity" of Table 2, with the corresponding entries: line 1 for GLASSO-based sparsity; line 2 for D-trace SCAD-based sparsity; line 3 for D-trace MCP-based sparsity; line 4 for D-trace LASSO-based sparsity; line 5 for D-trace aLASSO-based sparsity.
  - the vector variable "vMSFEa" (line 256): provides column 4 "MSFE" of Table 2, with the corresponding entries: line 1 for GLASSO-based MSFE; line 2 for D-trace SCAD-based MSFE; line 3 for D-trace MCP-based MSFE; line 4 for D-trace LASSO-based MSFE; line 5 for D-trace aLASSO-based MSFE; line 6 for VAR OLS-based MSFE; line 7 for LASSO VAR OLS-based MSFE; line 8 for Ridge VAR OLS-based MSFE.
  - the matrix variable "MCS" (line 273): provides column 5 "MCS" of Table 2 (MCS procedure based on the range statistics provided by Hansen, Lunde and Nason, 2011), where the indices displayed in the first column of "MCS" correspond to: 1 for GLASSO; 2 for SCAD D-trace; 3 for MCP D-trace; 4 for LASSO D-trace; 5 for aLASSO D-trace; 6 for OLS VAR; 7 for LASSO OLS VAR; 8 for Ridge OLS VAR.

**C. Description of the main functions**

**sparse_precision.m**:

Purpose of the function: estimate a sparsity-based precision matrix based on a matrix of observations, where the optimal tuning parameter is selected by cross-validation (procedure adapted to a time series context).
<p align="center">
[theta_est,lambda_opt]= sparse_precision(X,loss,lambda,method,a_scad,b_mcp)
</p>

Inputs:
- X: vector of observation.
- loss: 'Gaussian' or 'DTrace'; if 'Gaussian', the GLASSO algorithm is performed (only with LASSO).
- lambda: tuning parameter (grid of values).
- method (for D-trace loss): 'scad', 'mcp', 'alasso', 'lasso'.
- a_scad: value of the scad parameter.
- b_mcp: value of the mcp parameter.
   
Outputs:
- theta_est: vech(Theta), i.e. column vector that stacks the columns of the lower triangular part of Theta.
- lambda_opt: optimal tuning parameter value chosen by cross-validation; if no-cross validation is performed (i.e., the input lambda is a scalar), then lambda_opt = lambda.
 
**dtrace_admm.m**:

Purpose of the function: estimate a sparsity-based precision matrix with D-trace loss. 
 
 <p align="center">
 Theta = dtrace_admm(S,W,lambda,method,a_scad,b_mcp)
 </p>
 
The procedure is based on the ADMM algorithm proposed in Zhang and Zou (2014): the main code corresponds to Algorithm 2 of Zhang and Zou (2014); if the estimated Theta is not positive definite, then Algorithm 1 of Zhang and Zou (2014) is performed.For the SCAD and MCP penalty functions, the local linear approximation (LLA) algorithm of Zou and Li (2008) is implemented: the standard LASSO is run first; then the entries of this estimator enter the weight in the LLA algorithm: the weight is the first order derivative of the SCAD/MCP evaluated at the LASSO estimator.

Inputs:
- S: sample variance covariance matrix of the observations.
- W: matrix of weights active for the local linear approximation algorithm for SCAD and MCP; for SCAD and MCP, W is the LASSO estimator; for the LASSO penalization, W = $I_p$ the identity matrix.
- lambda: tuning parameter.
- method: 'scad', 'mcp', 'alasso', 'lasso'.
- a_scad: value of the scad parameter.
- b_mcp: value of the mcp parameter.

Outputs:
- Theta: sparse precision matrix estimated by the D-trace loss.
     
**graphicalLasso.m**:

Purpose of the function: estimate a sparsity-based precision matrix with Gaussian LASSO.
<p align="center">
[theta,W] = graphicalLasso(X,lambda,maxIt,tol)
</p>

The procedure performs the graphical LASSO proposed by Friedman et al. (2008). The LASSO is solved using the shooting algorithm (alternative LASSO methods can be considered).

Inputs:
- X: vector of observation. 
- lambda: tuning parameter.
- maxIt: maximum number of iterations (optional).
- tol: convergence tolerance level (optional).

Outputs:
- theta: vech(Theta), i.e. column vector that stacks the columns of the lower triangular part of Theta with Theta the sparse precision matrix estimated by the LASSO Gaussian loss.
- W: regularized covariance matrix estimator, W = $\text{Theta}^{-1}$.

**simulate_sparse_SVAR.m**:

Purpose of the function: generate a sparse $\Theta$ precision matrix, whose structure is deduced from sparse gVAR coefficients (see Section S.2 of the Online Supplement for specific details).
<p align="center">
[Theta,mAs,mSig_u,mB,mB0] = simulate_sparse_SVAR(N,P,sparsity,sparsity2,density,a1,a2)
</p>

Input:
- N: dimension of the vector of observations.
- P: number of desired VAR/SVAR lags.
- sparsity: desired sparsity degree in Theta_true, i.e., the number corresponding to the proportion of zero coefficients located in vech(Theta_true), the column vector that stacks the columns of the lower triangular part of Theta_true, excluding the diagonal coefficients. The proportion of zero entries is simply defined as sparsity/(q*(q-1)/2),  with q = N(P+1).
- sparsity2: controls for the minimum sparsity degree of the zeros in $B_0, B1, \ldots, B_p$: it refers to $\tilde{s}$ defined step (v) in the simulation procedure described in S.2 of the Online Supplement.
- density: controls the proportion of non-zero elements in the $q \times q$ matrix, where each entry is the sum of one or more normally distributed random samples.
- a1 and a2: coefficients controlling for the generating $C^{(1)}$ in step (ii) of the simulation procedure, i.e. $C^{(1)} = K + \beta I_q$, with $\beta \in \mathcal{U}([a1,a2])$ to ensure the positive-definiteness of $C^{(0)}, q = N(P+1)$.

Outputs: 
- Theta: $q \times q$ matrix precision matrix deduced from the gVAR model with q = N(P+1).
- mAs: VAR autoregressive coefficients deduced from Theta.
- mSig_u: variance-covariance matrix of the error system vector, deduced from $B_0$.
- mB: sparse SVAR autoregressive matrix coefficients $B_1, B_2,\ldots,B_p$.
- mB0: sparse SVAR matrix $B_0$.

**VAR_penalized.m**:

Purpose of the function: estimate the VAR process by LASSO/Ridge-penalized OLS, equation-by-equation.
<p align="center">
[theta_est,lambda_opt]= VAR_penalized(Y,X,lambda,method)
</p>

Inputs:
- Y: response variable.
- X: covariates.
- lambda: tuning parameter.
- method: 'ridge' or 'lasso'. *Note*: If 'lasso' is selected, then 'no intercept' is selected to get the estimate in line 114 of *real_data_analysis.m*, when executing the function *sparse_VAR.m*: the version 12.0 or above of the Statistics and Machine Learning Toolbox should be installed, otherwise an error will occur.

Outputs: 
- theta_est: penalized solution.
- lambda_opt: optimal tuning parameter value chosen by cross-validation; if no-cross validation is performed (i.e. the input lambda is a scalar), then lambda_opt = lambda.

# References

- Friedman, J., T. Hastie and R. Tibshirani (2008). Sparse inverse covariance estimation with the graphical lasso. *Biostatistics*, 9, 432-441.

- Hansen, P. R., A. Lunde and J. M. Nason (2011). The model confidence set. *Econometrica*, 79, 453-497.

- Heber, G., A. Lunde, N. Shephard and K. Sheppard (2009). Oxford-Man Institute's realized library. *Library Version: 0.3. Oxford-Man Institute: University of Oxford.* https://realized.oxford-man.ox.ac.uk/data (accessed: 21 August 2022).

- Zhang, T. and H. Zou (2014). Sparse precision matrix estimation via lasso penalized D-trace loss. *Biometrika*, 101, 103-120.

- Zou, H. and R. Li (2008). One-step sparse estimates in nonconcave penalized likelihood models. *The Annals of Statistics*, 36, 1509-1533.

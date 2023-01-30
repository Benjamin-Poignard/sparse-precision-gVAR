% Description of the main files to reproduce the simulation and real data
% results

%% Computer and Matlab specs:

% The Matlab code was run on a Mac-OS Apple M1 Ultra with 20 cores and 128
% GB Memory

% Version of Matlab software on which the code was run: 
% 9.12.0.1975300 (R2022a) Update 3
% The following toolboxes should equip the license
% Statistics and Machine Learning Toolbox               Version 12.3
% Financial Toolbox                                     Version 6.3       
% Econometrics Toolbox                                  Version 6.0
% Parallel Computing Toolbox                            Version 7.6

% => the code is run with parallel computation for the cross-validation
% procedure (selection of lambda, the tuning parameter)

%% Simulation experiments (Replication of Table 1)

% File to run for the simulated experiment described in Section 4 of the
% paper: "simulations.m" ==> provides the results in Table 1

% - Case N = 10, p = 2: both sample sizes T = 500 and T = 2000 
% are included

% - Case N = 20, p = 1: both sample sizes T = 500 and T = 2000 
% are included

% For each batch, a true Theta, denoted by "Theta_true" is simulated,
% providing "mAs_true" and "mSig_u_true", the true VAR coefficients deduced
% from "Theta_true" and the true "B0, B1,....,Bp"
% ==> refer to "simulate_sparse_SVAR.m" to generate theses parameters

% The main function to estimate a sparse precison matrix is 
% sparse_precision.m, where the GLASSO, D-trace with SCAD, MCP, LASSO and
% adaptive LASSO are implemented; the penalisation parameter \lambda_T is
% obtained by an out-of-sample cross-validation procedure

% For each N, p, T case, the matrix "Results" is produced, with lines and
% columns:
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

% For the sake of transparency and replication, all the results based on 
% the simulated experiments are reported in the files:
% - "results_N10_P2.m" for N = 10, p = 2 and sample size T = 500, 2000
% - "results_N20_P1.m" for N = 20, p = 1 and sample size T = 500, 2000
% In those files, "Results" produces the results reported in Table 1 

%% Real data experiment (Replication of Table 2 and Figure 1)

% File to run for the real data experiment described in Section 5 of the
% paper: "real_data_analysis.m" ==> provides Figure 1 (named "figure2_rk") and Table 2

% The raw data are in "data_rk.csv" which are accessible from Oxford-Man
% Institutes realised library (Heber et al., 2009):
% https://www.oxford-man.ox.ac.uk/research/realized-library/
% Note that the server of Oxfor-Man may experience interruptions

% The data mD2 used in the real data analysis is created in the file 
% "create_data_rk.m" ==> produce "data_rk.mat", which contains mD1 and mD2
% mD2 data: daily realized return volatility of the following indices:

% #   Index                  Country                 
% 1   S&P 500                 USA              
% 2   FTSE 100                UK                      
% 3   Nikkei 225              Japan                    
% 4   DAX                     Germany
% 5   Russel 2000             USA
% 6   All Ordinaries          Australia                
% 7   CAC 40                  France    
% 8   Hang Seng               Hong Kong                          
% 9   KOSPI                   South Korea                          
% 10  AEX Index               The Netherlands                  
% 11  Swiss Market Index      Switzerland                        
% 12  IBEX 35                 Spain              
% 13  S&P CNX Nifty           India                  
% 14  IPC Mexico              Mexico                      
% 15  Bovespa Index           Brazil                
% 16  S&P/TSX                 Canada                 
% 17  Euro STOXX 50           50 largest companies in the eurozone
% 18  FTSE MIB                Italy

% Note: Euro STOXX 50 index is comprised of the 50 larges companies  
% in the following 12 eurozone countries: Austria, Belgium, Finland, 
% France, Germany, Greece, Ireland, Italy, Luxembourg, the Netherlands, 
% Portugal, Spain

% Note: the raw data file "data_rk.csv" contains the index STI (FT Straits 
% Times Index, Singapore), but the starting date of observations is
% different from the other indices, and, thus, is not used in the real data
% analysis

% data_rk.mat contains mD1 and mD2 ==> only mD2 is used and contains the
% aforementioned variables. The raw data in mD2 starts from 2000/01/03 to 
% 2017/12/04. The real data experiment carried out in Section 5 is based on
% the sub-sample 2009/11/23-2017/12/04

% How to obtain the results in Table 2 and Figure 1:
% run "real_data_analysis.m": rolling windown estimation of the gVAR 
%           with sparse precision matrix (with GLASSO and D-trace with
%           SCAD, MCP, LASSO and adaptive LASSO) and of the VAR with LASSO 
%           OLS and Ridge OLS. The the code computes the MFSE and MCS 
%           statistics based on the estimates saved in "mR_rk_all_for.mat" 

%  208 in "real_data_analysis.m": change the index i to get the
% matrix coefficients B0, B1,... Bp estimated by the desired method:
% i==1: GLASSO; i==2: D-trace SCAD; i==3: D-trace MCP; 
% i==4: D-trace LASSO; i==5: D-trace aLASSO;  

% The one-step-ahead forecasts are contained in mR_rk_all_for.mat, with:
% - mF_lasso_g: Gaussian with LASSO penalization precision matrix
% - mF_scad_dt: D-Trace with SCAD penalization precision matrix
% - mF_mcp_dt: D-Trace with MCP penalization precision matrix
% - mF_lasso_dt: D-Trace with LASSO penalization precision matrix
% - mF_alasso_dt: D-Trace with adaptive LASSO penalization precision matrix
% - mF_est_lasso: VAR with LASSO penalization (equation by equation)
% - mF_est_ridge: VAR with Ridge penalization (no sparsity is fostered with Ridge)

%% Additional comments: 

% On the tuning parameter lambda: 
% set as c*sqrt(log(dim)/T), with T the sample size and dim the number of
% parameters to estimate (dim is of order [N(p+1)]^2); c is a
% user-specified grid of values: set as 0.001, 0.011, etc; In the real data
% analysis, c is set as 2.5:0.01:4 for SCAD/MCP as, empirically, we found
% that SCAD/MCP do not tend to select small values of lambda, contrary to
% LASSO-based methods (GLASSO and D-trace LASSO)
% For the simulated experiments, the grid is set as (0.001:0.01:3) for all
% penalization methods

% On GLASSO:
% the input "maxIt", the maximum number of iterations is set as default as
% 1e5. 

% On D-trace:
% the maximum number of iterations is set as default as
% 1e5. 
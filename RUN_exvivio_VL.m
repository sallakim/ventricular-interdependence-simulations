clear all
close all

addpath ExVivo_3
addpath Core

printon = 1; 

%% Load pseudodata structure 

data = makedatastructure;

%% Global parameters

ODE_TOL = 1e-8; 
data.gpars.ODE_TOL = ODE_TOL; 

%% Set volume loading

a_eta_Vtot = [1:.05:2];
data.a_eta_Vtot = a_eta_Vtot; 

%% Load in parameters parameters 

[adjpars,~,~,data] = parameters(data); 

%% Solve model 

[outputs,rout,J] = model_sol_4gammaopt(adjpars,data); 

%% Plot figures 

plot_exvivo_VL(outputs,data)
    
% elapsed_time = toc


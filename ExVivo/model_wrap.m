%--------------------------------------------------------------------------
%Used to fix some parameters and let the others vary (INDMAP) before
%solving ODE
%--------------------------------------------------------------------------

function [J,sol,rout] = model_wrap(pars,data)

%{
    This function wraps around model_sol to modify only the parameters
    optimized. For this optimization, gamma is reassigned. 
    Inputs: 
    pars    - vector of parameters to optimize 
    data    - input data structure with data and global parameters 
    Outpus: 
    J       - cost functional
    sol     - solution output structure from model_sol_exvivo.m 
    rout    - residual vector 
%}


data.gamma_opt = exp(pars); 
[tpars,~,~,data] = parameters(data); 

[sol,rout,J] = model_sol_gammaopt(tpars,data);

% only want to load in the most significant parameters, reassign
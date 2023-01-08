%--------------------------------------------------------------------------
%Used to fix some parameters and let the others vary (INDMAP) before
%solving ODE
%--------------------------------------------------------------------------

function [J,sol,rout] = model_wrap(pars,data)

% ALLPARS = data.gpars.ALLPARS;
% INDMAP  = data.gpars.INDMAP;
% 
% tpars = ALLPARS;
% tpars(INDMAP') = pars;

data.gamma_opt = exp(pars); 
[tpars,~,~,data] = parameters(data); 

[sol,rout,J] = model_sol_4gammaopt(tpars,data);

% only want to load in the most significant parameters, reassign
function outputs_loading = volumeloading_wrap(pars,data,eta_Vtot)


%{ 
    This function wraps around the volume loading code to assign the
    scalar determining loading of the total circulating volume in the model
    called eta_Vtot. 
    Inputs: 
    optpars         - vector of optimized parameters 
    data            - input data structure with data and global parameters 
    eta_Vtot        - scalar value for volume loading 
    Outputs: 
    outputs_loading - output structure from model_sol.m 
%} 

data.eta_Vtot = eta_Vtot; 
outputs_loading = model_sol(pars,data);

end
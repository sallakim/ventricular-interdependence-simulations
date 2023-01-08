function outputs_loading = volumeloading_wrap(pars,data,eta_Vtot)

data.eta_Vtot = eta_Vtot; 
outputs_loading = model_sol(pars,data);

end
clear all
close all

addpath ExVivo_3
addpath Core

printon = 1; 

% Print output figures from model 
printoutfigs_on = 1; % 0 - off; 1 - on 

%% Load pseudodata structure 

data = makedatastructure;

%% Global parameters

ODE_TOL = 1e-8; 
data.gpars.ODE_TOL = ODE_TOL; 

%% Set volume loading

a_eta_Vtot = [1:.1:2];
data.a_eta_Vtot = a_eta_Vtot; 

eta_EDV_vec = [0.5,1,1.5,2]; 

for i = 1:length(eta_EDV_vec)

    % Scale EDV
    data.EDV_LV = eta_EDV_vec(i) * data.EDV_LV;
    data.EDV_RV = eta_EDV_vec(i) * data.EDV_RV; 
    
    % Load in parameters
    [adjpars,~,~,data] = parameters(data); 
    
    % Solve model
    [outputs,~,~] = model_sol_4gammaopt(adjpars,data); 

    % Undo EDV scaling
    data.EDV_LV = 1/eta_EDV_vec(i) * data.EDV_LV; 
    data.EDV_RV = 1/eta_EDV_vec(i) * data.EDV_RV; 

    EDV_LV_mat(:,i) = outputs.EDPVRs.a_EDV_LV(:); 
    EDV_RV_mat(:,i) = outputs.EDPVRs.a_EDV_RV(:); 

    EDP_LV_mat(:,i) = outputs.EDPVRs.a_EDP_LV(:); 
    EDP_RV_mat(:,i) = outputs.EDPVRs.a_EDP_RV(:); 

end

outputs_EDV.eta_EDV_vec = eta_EDV_vec; 
outputs_EDV.EDV_LV_mat  = EDV_LV_mat; 
outputs_EDV.EDV_RV_mat  = EDV_RV_mat; 
outputs_EDV.EDP_LV_mat  = EDP_LV_mat; 
outputs_EDV.EDP_RV_mat  = EDP_RV_mat; 

%% Plot figures 

data.printoutfigs_on = printoutfigs_on;

plot_exvivo_VL_EDV(outputs_EDV,data)
    
% elapsed_time = toc
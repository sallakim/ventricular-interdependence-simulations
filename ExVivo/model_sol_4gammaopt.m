function [outputs,rout,J] = model_sol_4gammaopt(adjpars,data)

%% Parameters  

% "undo" log from parameters.m
adjpars = exp(adjpars); 

%% Unpack data structure 

a_eta_Vtot = data.a_eta_Vtot; 


%% Solve the model 

for q = 1:length(a_eta_Vtot)

    % Assign volume loading vector
    eta_Vtot = a_eta_Vtot(q); 

    % Get initial conditions 
    init = initialconditions(adjpars,data,eta_Vtot);

    % Run model 

    [~,o] = model_exvivo([],init,adjpars,data);

    V_LV = init(8) * 1e6; 
    V_RV = init(11) * 1e6; 
    P_LV = o(1,:) * 7.5;
    P_RV = o(2,:) * 7.5;

    % Arrays
    a_EDP_LV(q) = P_LV;
    a_EDV_LV(q) = V_LV;
    a_EDP_RV(q) = P_RV;
    a_EDV_RV(q) = V_RV;
end

%% Klotz curve calculation 

% Klotz et al. EDPVR Calculation  
An = 28; 
Bn = 3; 

% LV
EDP_LV_data = 5; 
EDV_LV_data = 125; 
EDP_RV_data = 5/4; 
EDV_RV_data = 125; 

V_0_LV  = EDV_LV_data * (0.6 - 0.006 * EDP_LV_data); 
V_30_LV = V_0_LV + (EDV_LV_data - V_0_LV) / ((EDP_LV_data / An)^(1/Bn));
Beta_LV = log(EDP_LV_data/30) / log(EDV_LV_data / V_30_LV); 
Alpha_LV = 30 / V_30_LV^Beta_LV; 
P_LV_EDPVR = Alpha_LV * a_EDV_LV.^Beta_LV;

V_0_RV  = EDV_RV_data * (0.6 - 0.006 * EDP_RV_data); 
V_30_RV = V_0_RV + (EDV_RV_data - V_0_RV) / ((EDP_RV_data / An)^(1/Bn));
Beta_RV = log(EDP_RV_data/30) / log(EDV_RV_data / V_30_RV); 
Alpha_RV = 30 / V_30_RV^Beta_RV; 
P_RV_EDPVR = Alpha_RV * a_EDV_RV.^Beta_RV;

EDPVRs.V_LV_EDPVR = a_EDV_LV;
EDPVRs.V_RV_EDPVR = a_EDV_RV; 
EDPVRs.P_LV_EDPVR = P_LV_EDPVR;  
EDPVRs.P_RV_EDPVR = P_RV_EDPVR; 

EDPVRs.a_EDV_LV = a_EDV_LV; 
EDPVRs.a_EDV_RV = a_EDV_RV; 
EDPVRs.a_EDP_LV = a_EDP_LV; 
EDPVRs.a_EDP_RV = a_EDP_RV; 

%% Outputs 

r1 = (a_EDP_LV - P_LV_EDPVR)./P_LV_EDPVR; 
r2 = (a_EDP_RV - P_RV_EDPVR)./P_RV_EDPVR; 
rout = [r1'; r2']; 

EDPVRs.rout = rout;
outputs.EDPVRs = EDPVRs; 

J = rout'*rout;

function [P_LV_EDPVR,P_RV_EDPVR,V] = makeKlotzcurve(EDV,EDP,V_EDPVR)

%{ 
    This function makes the single-beat approximated end-diastolic pressure
    -volume relationship from Klotz et al. (Am J Physiol Heart Circ Physiol
    (2006).) 

    Inputs 
    EDV         - vector of end-diastolic volumes for the LV and RV 
    EDP         - vector of end-systolic volumes for the LV and RV 
    V_EDPVR     - vector of discretized volumes as an input for the EDPVR
    
    Outputs: 
    P_LV_EDPVR  - vector of corresponding LV EDPVR pressures 
    P_RV_EDPVR  - vector of corresponding RV EDPVR pressures 

%} 

%% Parameters set from Klotz et al. 
An = 28; 
Bn = 3; 

%% LV 

EDP_LV_norm = EDP(1); 
EDV_LV_norm = EDV(1); 
    
V_0_LV  = EDV_LV_norm * (0.6 - 0.006 * EDP_LV_norm); 
V_30_LV = V_0_LV + (EDV_LV_norm - V_0_LV) / ((EDP_LV_norm / An)^(1/Bn));
Beta_LV = log(EDP_LV_norm/30) / log(EDV_LV_norm / V_30_LV); 
Alpha_LV = 30 / V_30_LV^Beta_LV; 

%% RV
    
EDP_RV_normal = EDP(2); 
EDV_RV_normal = EDV(2); 
    
V_0_RV  = EDV_RV_normal * (0.6 - 0.006 * EDP_RV_normal); 
V_30_RV = V_0_RV + (EDV_RV_normal - V_0_RV) / ((EDP_RV_normal / An)^(1/Bn));
Beta_RV = log(EDP_RV_normal/30) / log(EDV_RV_normal / V_30_RV); 
Alpha_RV = 30 / V_30_RV^Beta_RV; 
    
%% Outputs

P_LV_EDPVR = Alpha_LV * V_EDPVR.^Beta_LV;
P_RV_EDPVR = Alpha_RV * V_EDPVR.^Beta_RV;

V = [V_0_LV; V_30_LV; V_0_RV; V_30_RV]; 

end
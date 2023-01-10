function [data] = makedatastructure

%{
    This function compiles all of the data into a structure.
    Inputs: 
    data    - input data structure with data and global parameters 
    Outputs: 
    data    - reassigned data structure with data from this script included
%} 

% Heart rate (H) and heart period (T)
HR = 60;  
T  = 60/HR;  

% Set tspan to solve over the first 20 beats 
dt    = 0.001; 
tspan = 0:dt:20*T;

%% Load pseudodata

% Mean systolic and diastolic blood pressure 
SPbar = 120; %140
DPbar = 80; %90

% Total blood volume (mL)
Vtot = 5000; 

% Cardiac Output (m^3 s^(-1))
CO_LV = Vtot / 60;  

SV = CO_LV * T; 

% End-diastolic and end-systolic pressures and volumes 
EDP_LV = 5;         
EDP_RV = EDP_LV/4;           %EDP_lv * .25;  % Assume EDP_rv is 25% of EDP_lv 
EDV_LV = 125;      
EDV_RV = EDV_LV;        % Assume EDV_rv is equal to EDV_lv 
ESP_LV = 1.05*SPbar;    %125;       
ESP_RV = 25; 
ESV_LV = EDV_LV - SV; %50;       
ESV_RV = ESV_LV;        % Assume EDV_rv is equal to EDV_lv 

% fraction of unstressed volume in the compartments 
bvd_SA = .7; 
bvd_PA = .4;
bvd_SV = .9;
bvd_PV = .9;

% Load data into a structure and convert units to SI 
data.SPbar = SPbar / 7.5; 
data.DPbar = DPbar / 7.5; 
data.HR     = HR; 
data.T      = T; 
data.Vtot   = Vtot  * 1e-6; 
data.CO     = CO_LV * 1e-6; 
data.tspan  = tspan; 
data.dt     = dt;
data.EDP_LV = EDP_LV / 7.5; 
data.EDV_LV = EDV_LV * 1e-6;
data.ESP_LV = ESP_LV / 7.5; 
data.ESV_LV = ESV_LV * 1e-6; 
data.EDP_RV = EDP_RV / 7.5; 
data.EDV_RV = EDV_RV * 1e-6; 
data.ESP_RV = ESP_RV / 7.5; 
data.ESV_RV = ESV_RV * 1e-6; 

data.bvd_SA = bvd_SA; 
data.bvd_PA = bvd_PA; 
data.bvd_SV = bvd_SV; 
data.bvd_PV = bvd_PV; 

%% Experimental scaling factors

eta_Vtot     = 1;  % 2.5 for PSD
eta_k_pas_LV = 1; 
eta_k_pas_RV = 1; 
eta_k_act_LV = 1;  % 1.5 for PSD
eta_k_act_RV = 1; 

data.eta_Vtot = eta_Vtot;
data.eta_k_pas_LV = eta_k_pas_LV; 
data.eta_k_pas_RV = eta_k_pas_RV; 
data.eta_k_act_LV = eta_k_act_LV; 
data.eta_k_act_RV = eta_k_act_RV; 

end
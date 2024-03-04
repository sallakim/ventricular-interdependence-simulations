function [adjpars,UB,LB,data] = parameters(data)

 %{ 
    Assignment and/or nominal calculation of all model parameters. 
    
    Inputs: 
    data        - input data structure with data and global parameters 
    Outputs: 
    adjpars     - vector of adjustable parameters 
    UB          - vector of parameter upper bounds 
    LB          - vector of parameter lower bounds 
    data        - output data structure with new field assignments 
 %} 

%% Load in values from data structure 
% Blood pressures 
SPbar = data.SPbar; 
DPbar = data.DPbar; 

% Total blood volume (m^3) 
Vtot  = data.Vtot; 

% Cardiac output (m^3 s^(-1))
CO    = data.CO;

% End-diastolic and end-systolic pressures (kPa) and volumes (m^3) 
ESP_LV = data.ESP_LV; 
ESV_LV = data.ESV_LV; 
EDP_LV = data.EDP_LV; 
EDV_LV = data.EDV_LV; 

ESP_RV = data.ESP_RV; 
EDP_RV = data.EDP_RV;   
EDV_RV = data.EDV_RV;   

% Experimental loading factors
eta_k_pas_LV = data.eta_k_pas_LV; 
eta_k_pas_RV = data.eta_k_pas_RV; 
eta_k_act_LV = data.eta_k_act_LV; 
eta_k_act_RV = data.eta_k_act_RV; 

%% Pericardium 

Vh0 = 1.25*(EDV_LV + EDV_RV); 
s = 10; % increase / decrease to change pericardial pressures 

%% Blood volume distribution 

bvd_SA = data.bvd_SA;   bvd_PA = data.bvd_PA; 
bvd_SV = data.bvd_SV;   bvd_PV = data.bvd_PV; 

% sum total = 1.0, total fraction of blood in each compartment 
% d_LV = .03;             d_RV = .03;
% d_SA = .2;              d_PA = .07; 
% d_SV = .45;             d_PV = .20;

d_LV = .025;           d_RV = .025;
d_SA = .2;             d_PA = .05; 
d_SV = .45;             d_PV = .25;

d_sum_LV = d_PV + d_LV;
d_sum_RV = d_SV + d_RV;

d_LV = EDV_LV / Vtot; 
d_RV = EDV_RV / Vtot;
d_PV = d_sum_LV - d_LV;
d_SV = d_sum_RV - d_RV;

%% Calculate volumes (m^3)

% Total chamber volumes: bld volume distribution fraction x total volume 
V_LV_0 = d_LV*Vtot;      V_RV_0 = d_RV*Vtot; 
V_SA_0 = d_SA*Vtot;      V_PA_0 = d_PA*Vtot; 
V_SV_0 = d_SV*Vtot;      V_PV_0 = d_PV*Vtot;

% Unstressed volumes
V_SA_u = V_SA_0*bvd_SA;   V_PA_u = V_PA_0*bvd_PA; 
V_SV_u = V_SV_0*bvd_SV;   V_PV_u = V_PV_0*bvd_PV; 

% Stressed volumes
V_SA_s = V_SA_0 - V_SA_u;  V_PA_s = V_PA_0 - V_PA_u;  
V_SV_s = V_SV_0 - V_SV_u;  V_PV_s = V_PV_0 - V_PV_u; 

%% Pressures (kPa)
% Max/min pressure ratios from Boron book (mmHg converted to kPa)
% LV                    % SA                    % SV
P_LV_M  = 121 / 7.5;     P_SA_M   = 120 / 7.5;    P_SV_M   = 8 / 7.5;
                         P_SA_bar = 95  / 7.5;    P_SV_bar = 3 / 7.5; 
P_LV_m  = 1   / 7.5;     P_SA_m   = 80  / 7.5;    P_SV_m   = 2 / 7.5;
                                            
% RV                    % PA                    % PV
P_RV_M = 25 / 7.5;       P_PA_M   = 25 / 7.5;     P_PV_M   = 8  / 7.5; 
                         P_PA_bar = 15 / 7.5;     P_PV_bar = 3 / 7.5; 
P_RV_m = 2.5  / 7.5;     P_PA_m   = 10 / 7.5;    P_PV_m   = 2  / 7.5; 

% Scale pressures to SA pressures 
q_LV_M = P_LV_M/P_SA_M;    q_LV_m = P_LV_m/P_SA_m; 
q_RV_M = P_RV_M/P_SA_M;    q_RV_m = P_RV_m/P_SA_m; 
q_SA_M = P_SA_M/P_SA_M;    q_SA_m = P_SA_m/P_SA_m; 
q_SV_M = P_SV_M/P_SA_M;    q_SV_m = P_SV_m/P_SA_m; 
q_PA_M = P_PA_M/P_SA_M;    q_PA_m = P_PA_m/P_SA_m; 
q_PV_M = P_PV_M/P_SA_M;    q_PV_m = P_PV_m/P_SA_m; 

% Recapitulate nominal pressures for each chamber given SA pressure data    
P_LV_M = q_LV_M*SPbar;    P_LV_m = q_LV_m*DPbar;       
P_RV_M = q_RV_M*SPbar;    P_RV_m = q_RV_m*DPbar;       
P_SA_M = q_SA_M*SPbar;    P_SA_m = q_SA_m*DPbar;       
P_SV_M = q_SV_M*SPbar;    P_SV_m = q_SV_m*DPbar;
P_PA_M = q_PA_M*SPbar;    P_PA_m = q_PA_m*DPbar;       
P_PV_M = q_PV_M*SPbar;    P_PV_m = q_PV_m*DPbar;
                                            
%% Calculate elastances (kPa m^(-3)) and compliances (m^3 kPa^(-1))

% LV and RV minimal elastance parameters 
E_LV_m = P_LV_m / V_LV_0; 
E_RV_m = P_RV_m / V_RV_0; 

% Compliances
C_SA = V_SA_s/P_SA_M;
C_SV = V_SV_s/P_SV_M; 
C_PA = V_PA_s/P_PA_M; 
C_PV = V_PV_s/P_PV_M; 

%% Calculate resistances (kPa s m^(-3)) 
R_SA = (P_SA_M - P_SV_bar)/CO;
R_PA = (P_PA_M - P_PV_bar)/CO; 

% Transmural resistances 
R_tSA = 0.08 / 7.5e-6; % orig 0.05
R_tPA = 0.02 / 7.5e-6; % orginal tPA = tSA

% Valves (vlv)
R_m = (P_PV_bar - P_LV_m) / CO; %6e-2 / 7.5e-6; 
R_a = 5e-4 / 7.5e-6; 
R_t = (P_SV_bar - P_RV_m) / CO; %8e-3 / 7.5e-6; 
R_p = 5e-4 / 7.5e-6;

% if eta_k_pas_LV ~= 1 || eta_k_pas_RV ~= 1 || eta_k_act_LV ~= 1 || eta_k_act_RV ~= 1 
%     R_m = R_m/2;
%     R_t = R_t/2; 
% end 

%% Heart model parameters 

% Sarcomere length parameters (convert from µm to m)
Lsref   = 2    * 1e-6; 
Lsc0    = 1.51 * 1e-6; 
Lse_iso = 0.04 * 1e-6; 

% Sarcomere length shortening velocity (convert from µm s^(-1) to m s^(-1))
v_max   = .5*7 * 1e-6;    

% Passive stress steepness parameter  
gamma = 7.5; %7.38; 
if isfield(data,'opt_gamma')
    gamma = data.opt_gamma; 
end 

% Unit conversion factor 
nu_SL = 1e6; 

%% Calculate subject-specific reference midwall surface area (Amref) for LV, SEP, and RV

% Ventricular inner chamber radius 
r_LV_and_SEP = (EDV_LV * 3 / (4* pi))^(1/3); % use eqn of sphere to get radius
r_RV         = (EDV_RV * 3 / (4* pi))^(1/3); 

% Ventricle midwall radius (chamber radius (r) + 1/2 wall thickness (h)) where h is
h_LV_and_SEP = 8 * 1e-3; 
h_RV         = 4 * 1e-3; % 3-5 mm https://heart.bmj.com/content/early/2022/05/30/heartjnl-2021-320526
% Midwall radius 
r_m_LV_and_SEP = r_LV_and_SEP + h_LV_and_SEP/2; % m 
r_m_RV         = r_RV + h_RV/2; %m 

% Outer radius 
r_o_LV_and_SEP = r_LV_and_SEP + h_LV_and_SEP; % m 
r_o_RV         = r_RV + h_RV; %m 

% Midwall reference surface area 
Amref_LV_and_SEP = 4 * pi * (r_m_LV_and_SEP)^2; % assume sphere 
Am_RV            = 4 * pi * (r_m_RV)^2;

Amref_LV  = Amref_LV_and_SEP * 2/3; % Assume LV is 2/3 of LV+SEP 
Amref_SEP = Amref_LV_and_SEP * 1/3; % Assume SEP is 1/3 of LV+SEP
Amref_RV  = Am_RV;

%% Calculate patient-specific wall volume (Vw) for LV, SEP, and RV 

% Outer Ventricle volume 
Vw_chamber_LV_and_SEP = 4/3 * pi * r_o_LV_and_SEP^3;  
Vw_chamber_RV         = 4/3 * pi * r_o_RV^3; 

% Ventricular wall volume 
Vw_LV_and_SEP = Vw_chamber_LV_and_SEP - EDV_LV; 
Vw_RV         = Vw_chamber_RV - EDV_RV;  

Vw_LV  = Vw_LV_and_SEP * 2/3; % Assume LV is 2/3 of LV+SEP 
Vw_SEP = Vw_LV_and_SEP * 1/3; % Assume SEP is 1/3 of LV+SEP 

%% Approximations for initial displacements and Amref_rv in end-diastole 

% Initialize diastolic displacement values 
xm_LV_d_0  = 5 * 1e-2; 
xm_SEP_d_0 = 2 * 1e-2; 
xm_RV_d_0  = 6 * 1e-2; 
ym_d_0     = 3 * 1e-2; 

x0 = [xm_LV_d_0; 
    xm_SEP_d_0; 
    xm_RV_d_0;
    ym_d_0; 
    Amref_RV; 
    ]; 
x0 = log(x0); % log-scale the initial values 

% Inputs for calculating displacements 
Vw    = [Vw_LV,Vw_SEP,Vw_RV]; 
Amref = [Amref_LV,Amref_SEP]; 

% Assume end-systolic sarcomere length is 1.6 um 
SL_d    = 2; %um 

opts = optimoptions('fsolve','Display','none',...
    'MaxFunctionEvaluations',2e3,'Algorithm','levenberg-marquardt'); 
[fnew0,~] = fsolve(@(x) calc_xm_ym(x,Lsref,Vw,Amref,SL_d,EDV_LV,0),x0,opts); 
fnew0 = exp(fnew0);

% Outputs / Diastolic displacements (m)  
xm_LV_d  = fnew0(1);
xm_SEP_d = fnew0(2);
xm_RV_d  = fnew0(3);
ym_d     = fnew0(4);
Amref_RV = fnew0(5); 

% Set initial conditions for displacements to be used in the
% initialconditions.m script 
deformation.xm_LV_0  = xm_LV_d; 
deformation.xm_SEP_0 = xm_SEP_d; 
deformation.xm_RV_0  = xm_RV_d; 
deformation.ym_0     = ym_d; 

data.deformation = deformation; 

%% Calculate passive stress parameters (k_pas) for LV and RV in end-diastole  

% Midwall surface area (m^2)
Am_LV_d  = pi * (xm_LV_d^2  + ym_d^2);
Am_RV_d  = pi * (xm_RV_d^2  + ym_d^2);

% Midwall curvature (m^(-1))
Cm_LV_d  = -2 * xm_LV_d  / (xm_LV_d^2  + ym_d^2);
Cm_RV_d  = -2 * xm_RV_d  / (xm_RV_d^2  + ym_d^2);

% Midwall ratio (dimensionless) 
z_LV_d   = 3 * Cm_LV_d  * Vw_LV  / (2 * Am_LV_d); 
z_RV_d   = 3 * Cm_RV_d  * Vw_RV  / (2 * Am_RV_d); 

% Instantaneous sarcomere length (um) in end-diastole
Ls_LV_d  = SL_d / nu_SL; 
Ls_RV_d  = SL_d / nu_SL;  

% Passive stress 
sigma_pas_LV_d = (nu_SL*(Ls_LV_d - Lsc0))^gamma; 
sigma_pas_RV_d = (nu_SL*(Ls_RV_d - Lsc0))^gamma; 

% Dimensionless combination function
Gamma_LV_d = -(2 / 3) * z_LV_d * (1 + (1 / 3) * z_LV_d^2 + (1 / 5) * z_LV_d^4);
Gamma_RV_d = -(2 / 3) * z_RV_d * (1 + (1 / 3) * z_RV_d^2 + (1 / 5) * z_RV_d^4);

% Passive stress scaling parameters (kPa)
k_pas_LV = eta_k_pas_LV*EDP_LV / (Gamma_LV_d * sigma_pas_LV_d);
k_pas_RV = eta_k_pas_RV*EDP_RV / (Gamma_RV_d * sigma_pas_RV_d);


%% Approximations for initial displacements and Amref_rv in end-systole 

% Initialize systolic displacements values 
xm_LV_d_0  = 5 * 1e-2; 
xm_SEP_d_0 = 2 * 1e-2; 
xm_RV_d_0  = 6 * 1e-2; 
ym_d_0     = 3 * 1e-2; 

x0 = [xm_LV_d_0; 
    xm_SEP_d_0; 
    xm_RV_d_0;
    ym_d_0; 
    ]; 
x0 = log(x0); % log-scale the initial values 

Amref = [Amref_LV,Amref_SEP, Amref_RV]; 

opts = optimoptions('fsolve','Display','none',...
    'MaxFunctionEvaluations',2e3,'Algorithm','levenberg-marquardt'); 
[fnew1,~] = fsolve(@(x) calc_xm_ym(x,Lsref,Vw,Amref,[],ESV_LV,1),x0,opts); 
fnew1 = exp(fnew1);

% Outputs / Systolic displacements (m)
xm_LV_s  = fnew1(1);
xm_RV_s  = fnew1(3);
ym_s     = fnew1(4);

%% Calculate active stress parameters (k_act) for LV and RV in end-systole 

% Midwall surface area (m^2)
Am_LV_s  = pi * (xm_LV_s^2  + ym_s^2);
Am_RV_s  = pi * (xm_RV_s^2  + ym_s^2);

% Midwall curvature (m^(-1))
Cm_LV_s  = - 2 * xm_LV_s  / (xm_LV_s^2  + ym_s^2);
Cm_RV_s  = - 2 * xm_RV_s  / (xm_RV_s^2  + ym_s^2);

% Midwall ratio (dimensionless)  
z_LV_s   = 3 * Cm_LV_s  * Vw_LV  / (2 * Am_LV_s); 
z_RV_s   = 3 * Cm_RV_s  * Vw_RV  / (2 * Am_RV_s); 

% Myofiber strain (dimensionless)
eps_LV_s  = 0.5 * log( Am_LV_s  / Amref_LV  ) - (1/12) * z_LV_s^2  - 0.019 * z_LV_s^4; 
eps_RV_s  = 0.5 * log( Am_RV_s  / Amref_RV  ) - (1/12) * z_RV_s^2  - 0.019 * z_RV_s^4; 

% Sarcomere length (m)
Ls_LV_s  = Lsref * exp(eps_LV_s); 
Ls_RV_s  = Lsref * exp(eps_RV_s); 

% Activation function 
y_v = .6; % set to 60% in end-systole 

% Active stress 
sigma_act_LV_s = y_v * (Ls_LV_s  - Lsc0) / 1e-6; 
sigma_act_RV_s = y_v * (Ls_RV_s  - Lsc0) / 1e-6; 

% Dimensionless combination function 
Gamma_LV_s = - (2 / 3) * z_LV_s * (1 + (1 / 3) * z_LV_s^2 + (1 / 5) * z_LV_s^4);
Gamma_RV_s = - (2 / 3) * z_RV_s * (1 + (1 / 3) * z_RV_s^2 + (1 / 5) * z_RV_s^4);

% Active stress scaling parameters (kPa) 
k_act_LV = eta_k_act_LV*ESP_LV / (Gamma_LV_s * sigma_act_LV_s);
k_act_RV = eta_k_act_RV*ESP_RV / (Gamma_RV_s * sigma_act_RV_s);

%% Time-varying elastance model parameters 

% Percentage of cardiac cycle 
k_TS = 0.3; % Beginning of cardiac cycle to maximal systole  
k_TR = 0.3; 

%% Outputs

adjpars = [C_SA; C_SV; C_PA; C_PV;                      % 1-4
    R_SA; R_tSA; R_PA; R_tPA;                           % 5-8
    R_m; R_a; R_t; R_p;                 % 9-12
    Amref_LV; Amref_SEP; Amref_RV;                      % 13-15
    Vw_LV; Vw_SEP; Vw_RV;                               % 16-18
    k_pas_LV; k_pas_RV; k_act_LV; k_act_RV;             % 19-22
    gamma;                                              % 23                                          % 28
    ]; 

UB = adjpars*10; 
LB = adjpars/10;

adjpars = log(adjpars); 
UB = log(UB);
LB = log(LB);

fixpars = [P_LV_m; P_SA_m; P_SV_m; P_RV_m; P_PA_m; P_PV_m;          % 1-6
        Lsref; Lsc0; Lse_iso;                                       % 7-9
        v_max;                                                      % 10
        V_SA_u; V_SV_u; V_PA_u; V_PV_u;                             % 11-14
        k_TS; k_TR;                                                 % 15-16
        E_LV_m; E_RV_m;                                             % 17-18
        nu_SL;                                                      % 19
        Vh0; s;                                                     % 20-21
    ]; 
    
data.fixpars = fixpars; 

end 
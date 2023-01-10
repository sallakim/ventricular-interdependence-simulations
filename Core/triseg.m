function f = triseg(x,pars,data,init)

%{ 
    This function contains the equations to calculate consistent initial
    conditions for the DAE. 
    Inputs: 
    x           - vector of states 
    pars        - vector of adjustable parameters 
    data        - input data structure with data and global parameters 
    init        - vector of initial conditions 
    Outputs: 
    f           - vector of equations for the root finder 
%} 

x = exp(x); 

HR = data.HR; 

fixpars = data.fixpars;

%% Adjustable Parameters

Amref_LV  = pars(13); 
Amref_SEP = pars(14); 
Amref_RV  = pars(15); 

% Midwall volume (m^3)
Vw_LV  = pars(16); 
Vw_SEP = pars(17); 
Vw_RV  = pars(18); 

% Force scaling factors (kPa) 
k_pas_LV = pars(19);
k_pas_RV = pars(20);
k_act_LV = pars(21); 
k_act_RV = pars(22); 

% Passive stress steepness parameter (dimensionless) 
gamma = pars(23); 

%% Fixed parameters

% Sarcomere length parameters (m)
Lsref   = fixpars(7);
Lsc0    = fixpars(8); 
Lse_iso = fixpars(9); 

% Percentage of cardiac cycle 
k_TS = fixpars(15); % Beginning of cardiac cycle to maximal systole  
k_TR = fixpars(16); % Relaxation (maximal systole to baseline) 

% Unit conversion factor 
nu_SL = fixpars(19); 

%% Variables 

xm_LV  = x(1); 
xm_SEP = x(2); 
xm_RV  = x(3); 
ym     = x(4); 

% Contractile element length 
Lsc_LV  = x(5); 
Lsc_SEP = x(6); 
Lsc_RV  = x(7); 

% Volumes 
V_LV = init(8); 
V_RV = init(11);

%% Activation function 

T = 60/HR; 
TS = k_TS * T; 
TR = k_TR * T; 

t = 0; 
tc_v = mod(t,T);
if tc_v >= 0 && tc_v < TS 
    y_v = 0.5*(1 - cos(pi*tc_v/TS)); 
elseif tc_v >= TS && tc_v < TR + TS 
    y_v = 0.5*(1 + cos(pi*(tc_v - TS)/TR)); 
else
    y_v = 0; 
end  

%% Heart and sarcomere model 

% Volume of spherical cap formed by midwall surface (m^3)
Vm_LV  = -(pi / 6) * xm_LV  * (xm_LV^2  + 3 * ym^2); 
Vm_SEP =  (pi / 6) * xm_SEP * (xm_SEP^2 + 3 * ym^2); 
Vm_RV  =  (pi / 6) * xm_RV  * (xm_RV^2  + 3 * ym^2); 

% Surface area of midwall surface (m^2) 
Am_LV  = pi * (xm_LV^2  + ym^2);
Am_SEP = pi * (xm_SEP^2 + ym^2); 
Am_RV  = pi * (xm_RV^2  + ym^2); 

% Curvature of midwall surface (m^(-1))
Cm_LV  = -2 * xm_LV  / (xm_LV^2  + ym^2);
Cm_SEP =  2 * xm_SEP / (xm_SEP^2 + ym^2);
Cm_RV  =  2 * xm_RV  / (xm_RV^2  + ym^2);

% Ratio of wall thickness to midwall radius of curvature (dimensionless)
z_LV  = 3 * Cm_LV  * Vw_LV  / (2 * Am_LV); 
z_SEP = 3 * Cm_SEP * Vw_SEP / (2 * Am_SEP); 
z_RV  = 3 * Cm_RV  * Vw_RV  / (2 * Am_RV);

% Myofiber strain (dimensionless)
eps_LV  = 0.5 * log( Am_LV  / Amref_LV  ) - (1/12) * z_LV^2  - 0.019 * z_LV^4; 
eps_SEP = 0.5 * log( Am_SEP / Amref_SEP ) - (1/12) * z_SEP^2 - 0.019 * z_SEP^4; 
eps_RV  = 0.5 * log( Am_RV  / Amref_RV  ) - (1/12) * z_RV^2  - 0.019 * z_RV^4; 

% Sarcomere length (m)
Ls_LV  = Lsref * exp(eps_LV); 
Ls_SEP = Lsref * exp(eps_SEP); 
Ls_RV  = Lsref * exp(eps_RV); 

% Passive stress (kPa) 
sigma_pas_LV  =  k_pas_LV * ((Ls_LV - Lsc0)*nu_SL)^gamma; 
sigma_pas_SEP =  k_pas_LV * ((Ls_SEP - Lsc0)*nu_SL)^gamma; 
sigma_pas_RV  =  k_pas_RV * ((Ls_RV - Lsc0)*nu_SL)^gamma; 

% Active stress (kPa)
sigma_act_LV  = k_act_LV * y_v  * (Lsc_LV  - Lsc0) * nu_SL * (Ls_LV  - Lsc_LV)  / Lse_iso; 
sigma_act_SEP = k_act_LV * y_v  * (Lsc_SEP - Lsc0) * nu_SL * (Ls_SEP - Lsc_SEP) / Lse_iso;
sigma_act_RV  = k_act_RV * y_v  * (Lsc_RV  - Lsc0) * nu_SL * (Ls_RV  - Lsc_RV)  / Lse_iso;

% Total stress (kPa)
sigma_LV  = sigma_act_LV  + sigma_pas_LV; 
sigma_SEP = sigma_act_SEP + sigma_pas_SEP; 
sigma_RV  = sigma_act_RV  + sigma_pas_RV; 

% Representative midwall tension (mmHg cm^(-2))
Tm_LV  = (Vw_LV  * sigma_LV  / (2 * Am_LV))  * (1 + (z_LV^2)/3  + (z_LV^4)/5); 
Tm_SEP = (Vw_SEP * sigma_SEP / (2 * Am_SEP)) * (1 + (z_SEP^2)/3 + (z_SEP^4)/5); 
Tm_RV  = (Vw_RV  * sigma_RV  / (2 * Am_RV))  * (1 + (z_RV^2)/3  + (z_RV^4)/5);

% Axial midwall tension component 
Tx_LV  = - Tm_LV  * 2 * xm_LV  * ym / (xm_LV^2  + ym^2); 
Tx_SEP =   Tm_SEP * 2 * xm_SEP * ym / (xm_SEP^2 + ym^2); 
Tx_RV  =   Tm_RV  * 2 * xm_RV  * ym / (xm_RV^2  + ym^2); 

% Radial midwall tension component 
Ty_LV  = Tm_LV  * (-xm_LV^2  + ym^2) / (xm_LV^2  + ym^2); 
Ty_SEP = Tm_SEP * (-xm_SEP^2 + ym^2) / (xm_SEP^2 + ym^2); 
Ty_RV  = Tm_RV  * (-xm_RV^2  + ym^2) / (xm_RV^2  + ym^2);

%% System of equations 

f1 = -V_LV - 0.5 * Vw_LV - 0.5 * Vw_SEP + Vm_SEP - Vm_LV; 
f2 = Tx_LV + Tx_SEP + Tx_RV;
f3 = V_RV + 0.5 * Vw_RV + 0.5 * Vw_SEP + Vm_SEP - Vm_RV;
f4 = Ty_LV + Ty_SEP + Ty_RV; 

f5 = (Ls_LV  - Lsc_LV)  /Lse_iso - 1;
f6 = (Ls_SEP - Lsc_SEP) /Lse_iso - 1;
f7 = (Ls_RV  - Lsc_RV)  /Lse_iso - 1;

f = [f1; f2; f3; f4; f5; f6; f7]; 
end 


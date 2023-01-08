function [dxdt, outputs] = model(t,x,pars,data) 

fixpars = data.fixpars;

if isfield(data,'eta_EDV') == 1
    eta_EDV = data.eta_EDV; 
else
    eta_EDV = 1; 
end

%% Adjustable parameters

% Compliance (m^3 kPa^(-1))
C_SA = pars(1); 
C_SV = pars(2); 
C_PA = pars(3); 
C_PV = pars(4); 

% Resistance (kPa s m^(-3))
R_SA = pars(5); 
R_tSA = pars(6); 
R_PA = pars(7); 
R_tPA = pars(8); 

R_m_vlv = pars(9); 
R_a_vlv = pars(10); 
R_t_vlv = pars(11); 
R_p_vlv = pars(12); 

% Midwall reference surface area (m^2)
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

% sarcomere length shortening velocity
v_max = fixpars(10); % m s^(-1) 

% Percentage of cardiac cycle 
k_TS = fixpars(15); % Beginning of cardiac cycle to maximal systole  
k_TR = fixpars(16); % Relaxation (maximal systole to baseline) 

% Unit conversion factor 
nu_SL = fixpars(19); 

Vh0 = fixpars(20);
s = fixpars(21);

%% Variables 

% Axial distance of midwall junction (m)
xm_LV  = x(1); 
xm_SEP = x(2); 
xm_RV  = x(3);

% Radial distance of midwall junction (m)
ym = x(4); 

% Contractile element length (m)
Lsc_LV  = x(5); 
Lsc_SEP = x(6); 
Lsc_RV  = x(7); 

% Volumes (m^3) 
V_LV = x(8); 
V_SA = x(9); 
V_SV = x(10);
V_RV = x(11);
V_PA = x(12); 
V_PV = x(13); 

%% Activation function

y_lv = activation(t,data); 

%% Heart model

% Volume of spherical cap formed by midwall surface (m^3)
Vm_LV  = -(pi / 6) * xm_LV  * (xm_LV^2  + 3 * ym^2); 
Vm_SEP = (pi / 6) * xm_SEP * (xm_SEP^2 + 3 * ym^2); 
Vm_RV  = (pi / 6) * xm_RV  * (xm_RV^2  + 3 * ym^2); 

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
sigma_pas_LV  =  ((Ls_LV - Lsc0)*nu_SL)^gamma; 
sigma_pas_SEP =  ((Ls_SEP - Lsc0)*nu_SL)^gamma; 
sigma_pas_RV  =  ((Ls_RV - Lsc0)*nu_SL)^gamma; 

% Active stress (kPa)
sigma_act_LV  = y_lv  * (Lsc_LV  - Lsc0) * nu_SL * (Ls_LV  - Lsc_LV)  / Lse_iso; 
sigma_act_SEP = y_lv  * (Lsc_SEP - Lsc0) * nu_SL * (Ls_SEP - Lsc_SEP) / Lse_iso;
sigma_act_RV  = y_lv  * (Lsc_RV  - Lsc0) * nu_SL * (Ls_RV  - Lsc_RV)  / Lse_iso;

% Scaled stressed (kPa) 
sigma_pas_LV  = k_pas_LV * sigma_pas_LV; 
sigma_pas_SEP = k_pas_LV * sigma_pas_SEP; 
sigma_pas_RV  = k_pas_RV * sigma_pas_RV; 

sigma_act_LV  = k_act_LV * sigma_act_LV; 
sigma_act_SEP = k_act_LV * sigma_act_SEP; 
sigma_act_RV  = k_act_RV * sigma_act_RV; 

% Total stress (kPa)
sigma_LV  = sigma_act_LV  + sigma_pas_LV; 
sigma_SEP = sigma_act_SEP + sigma_pas_SEP; 
sigma_RV  = sigma_act_RV  + sigma_pas_RV; 

% Representative midwall tension (kPa m)
Tm_LV  = (Vw_LV  * sigma_LV  / (2 * Am_LV))  * (1 + (z_LV^2)/3  + (z_LV^4)/5); 
Tm_SEP = (Vw_SEP * sigma_SEP / (2 * Am_SEP)) * (1 + (z_SEP^2)/3 + (z_SEP^4)/5); 
Tm_RV  = (Vw_RV  * sigma_RV  / (2 * Am_RV))  * (1 + (z_RV^2)/3  + (z_RV^4)/5);

% Axial midwall tension component (kPa m)
Tx_LV  = - Tm_LV  * 2 * xm_LV  * ym / (xm_LV^2  + ym^2); 
Tx_SEP = Tm_SEP * 2 * xm_SEP * ym / (xm_SEP^2 + ym^2); 
Tx_RV  = Tm_RV  * 2 * xm_RV  * ym / (xm_RV^2  + ym^2); 

% Radial midwall tension component (kPa m)
Ty_LV  = Tm_LV  * (-xm_LV^2  + ym^2) / (xm_LV^2  + ym^2); 
Ty_SEp = Tm_SEP * (-xm_SEP^2 + ym^2) / (xm_SEP^2 + ym^2); 
Ty_RV  = Tm_RV  * (-xm_RV^2  + ym^2) / (xm_RV^2  + ym^2);

% Ventricular pressure (kPa)
ptrans_LV = 2 * Tx_LV / ym; 
ptrans_RV = 2 * Tx_RV / ym; 
P_LV      = -ptrans_LV; 
P_RV      = ptrans_RV; 

% Calculate Pericardial Pressure 
Vh = V_LV + V_RV; 
P_peri = 0; %exp(s*(Vh/Vh0 - 1))/7.5; % convert mmHg to kPa 
P_LV = P_peri + P_LV; 
P_RV = P_peri + P_RV; 

%% Lumped circulatory model 
% Pressure (kPa)
P_SV = V_SV / C_SV; 
P_PV = V_PV / C_PV;  

% When aortic valve is closed 
Q_a = 0; 
P_SA = (R_SA*V_SA + C_SA*P_SV*R_tSA + C_SA*Q_a*R_SA*R_tSA)/(C_SA*(R_SA + R_tSA)); 
Q_SA = (V_SA - C_SA*P_SV + C_SA*Q_a*R_tSA)/(C_SA*(R_SA + R_tSA)); 

% When aortic valve is open 
if (P_SA < P_LV) * (V_LV > 0)
    Q_a = -(R_SA*V_SA - C_SA*P_LV*R_SA - C_SA*P_LV*R_tSA + C_SA*P_SV*R_tSA)/(C_SA*(R_a_vlv*R_SA + R_a_vlv*R_tSA + R_SA*R_tSA)); 
    Q_SA = (R_a_vlv*V_SA - C_SA*P_SV*R_a_vlv + C_SA*P_LV*R_tSA - C_SA*P_SV*R_tSA)/(C_SA*(R_a_vlv*R_SA + R_a_vlv*R_tSA + R_SA*R_tSA)); 
    P_SA = (R_a_vlv*R_SA*V_SA + C_SA*P_SV*R_a_vlv*R_tSA + C_SA*P_LV*R_SA*R_tSA)/(C_SA*(R_a_vlv*R_SA + R_a_vlv*R_tSA + R_SA*R_tSA)); 
end
Q_a = max(Q_a,0); 

% When pulmonary valve is closed 
Q_p = 0; 
P_PA = (R_PA*V_PA + C_PA*P_PV*R_tPA + C_PA*Q_p*R_PA*R_tPA)/(C_PA*(R_PA + R_tPA)); 
Q_PA = (V_PA - C_PA*P_PV + C_PA*Q_p*R_tPA)/(C_PA*(R_PA + R_tPA)); 

% When pulmonary valve is open 
if (P_PA < P_RV) * (V_RV > 0) 
    Q_p = -(R_PA*V_PA - C_PA*P_RV*R_PA + C_PA*P_PV*R_tPA - C_PA*P_RV*R_tPA)/(C_PA*(R_PA*R_p_vlv + R_PA*R_tPA + R_p_vlv*R_tPA)); 
    Q_PA = (R_p_vlv*V_PA - C_PA*P_PV*R_p_vlv - C_PA*P_PV*R_tPA + C_PA*P_RV*R_tPA)/(C_PA*(R_PA*R_p_vlv + R_PA*R_tPA + R_p_vlv*R_tPA)); 
    P_PA = (R_PA*R_p_vlv*V_PA + C_PA*P_PV*R_p_vlv*R_tPA + C_PA*P_RV*R_PA*R_tPA)/(C_PA*(R_PA*R_p_vlv + R_PA*R_tPA + R_p_vlv*R_tPA));  
end 
Q_p = max(Q_p, 0); 

Q_m = max((P_PV - P_LV) / R_m_vlv, 0); 
Q_t = max((P_SV - P_RV) / R_t_vlv, 0); 

%% ODEs

% 1 - 4
dxm_LV  = -V_LV - 0.5 * Vw_LV - 0.5 * Vw_SEP + Vm_SEP - Vm_LV; 
dxm_SEP = Tx_LV + Tx_SEP + Tx_RV;
dxm_RV  = V_RV + 0.5 * Vw_RV + 0.5 * Vw_SEP + Vm_SEP - Vm_RV;
dym     = Ty_LV + Ty_SEp + Ty_RV; 

% 5 - 7
dLsc_LV  = ((Ls_LV  - Lsc_LV)  / Lse_iso - 1) * v_max;
dLsc_SEP = ((Ls_SEP - Lsc_SEP) / Lse_iso - 1) * v_max;
dLsc_RV  = ((Ls_RV  - Lsc_RV)  / Lse_iso - 1) * v_max;

% 8 - 14
dV_LV = Q_m - Q_a; 
dV_SA = Q_a - Q_SA; 
dV_SV = Q_SA    - Q_t; 
dV_RV = Q_t - Q_p; 
dV_PA = Q_p - Q_PA; 
dV_PV = Q_PA    - Q_m; 

dxdt = [dxm_LV; dxm_SEP; dxm_RV; dym;           % 1-4
    dLsc_LV; dLsc_SEP; dLsc_RV;                 % 5-7
    dV_LV; dV_SA; dV_SV; dV_RV; dV_PA; dV_PV;   % 8-13
    ]; 

outputs = [P_LV; P_SA; P_SV; P_RV; P_PA; P_PV;  % 1-6
    Vm_LV; Vm_SEP; Vm_RV;                       % 7-9
    Am_LV; Am_SEP; Am_RV;                       % 10-12
    Cm_LV; Cm_SEP; Cm_RV;                       % 13-15
    eps_LV; eps_SEP; eps_RV;                    % 16-18
    sigma_pas_LV; sigma_pas_SEP; sigma_pas_RV;  % 19-21
    sigma_act_LV; sigma_act_SEP; sigma_act_RV;  % 22-24
    sigma_LV; sigma_SEP; sigma_RV;              % 25-27
    Q_m; Q_a; Q_t; Q_p;         % 28-31
    Q_SA; Q_PA;                                 % 32-33
    Tm_LV; Tm_SEP; Tm_RV;                       % 34-36
    y_lv;                                        % 37
    P_peri;                                     % 38
    ];

end 
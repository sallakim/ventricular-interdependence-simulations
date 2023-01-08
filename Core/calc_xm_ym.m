function f = calc_xm_ym(x,Lsref,Vw,Amref,SL,V_LV,flag)
%{ 
This function is called in the parameters and calculates estimates for
the radial and axial midwall displacements for the LV, SEP, and RV. 
Inputs: 
    x       State vector 
    Vw      Vector of midwall volumes 
    Amref   Vector of reference surface areas 
    SL      Sarcomere length assumption (in diastole only)
    V_lv    LV volume during either end-diastole or end-systole 
    flag    0 - diastole; 1 - systole 
During diastole, this code estimates Amref_rv as a state and incorporates
an approximation for the diastolic sarcomere length. The diastolic
calculations are a fully determined system. 
During systole, this code uses the Amref_rv determined in diastole and 
solves an underdetermined system consisting of 2 equations.
%}

%% Parameters 

% Midwall volumes 
Vw_LV  = Vw(1); 
Vw_SEP = Vw(2); 
Vw_RV  = Vw(3); 

% Midwall references surface areas
Amref_LV  = Amref(1); 
Amref_SEP = Amref(2); 

% For end-systole  
if flag == 1 
    Amref_RV = Amref(3); 
end 

%% States

x = exp(x); 

xm_LV    = x(1); 
xm_SEP   = x(2); 
xm_RV    = x(3);
ym       = x(4); 

% For end-diastole
if flag == 0 
    Amref_RV = x(5); 
end 

%% Equations

% Volume of spherical cap formed by midwall surface (m^3)
Vm_LV  = -(pi / 6) * xm_LV  * (xm_LV^2  + 3 * ym^2); 
Vm_SEP = (pi / 6) * xm_SEP * (xm_SEP^2 + 3 * ym^2); 
Vm_RV  = (pi / 6) * xm_RV  * (xm_RV^2  + 3 * ym^2); 

% Midwall surface area (m^2)
Am_LV  = pi * (xm_LV^2  + ym^2);
Am_SEP = pi * (xm_SEP^2  + ym^2);
Am_RV  = pi * (xm_RV^2  + ym^2);

% Midwall curvature (m^(-1))
Cm_LV  = -2 * xm_LV  / (xm_LV^2  + ym^2);
Cm_SEP = -2 * xm_SEP / (xm_SEP^2  + ym^2);
Cm_RV  = -2 * xm_RV  / (xm_RV^2  + ym^2);

% Midwall ratio (dimensionless) 
z_LV   = 3 * Cm_LV  * Vw_LV  / (2 * Am_LV); 
z_SEP  = 3 * Cm_SEP * Vw_SEP / (2 * Am_SEP); 
z_RV   = 3 * Cm_RV  * Vw_RV  / (2 * Am_RV); 

% Strain (dimensionless) 
eps_LV  = 0.5 * log( Am_LV  / Amref_LV  ) - (1/12) * z_LV^2  - 0.019 * z_LV^4; 
eps_SEP = 0.5 * log( Am_SEP / Amref_SEP ) - (1/12) * z_SEP^2 - 0.019 * z_SEP^4; 
eps_RV  = 0.5 * log( Am_RV  / Amref_RV  ) - (1/12) * z_RV^2  - 0.019 * z_RV^4; 

% Instantaneous sarcomere length (um) 
Ls_LV  = Lsref * exp(eps_LV); 
Ls_SEP = Lsref * exp(eps_SEP); 
Ls_RV  = Lsref * exp(eps_RV); 

Ls_LV  = Ls_LV  * 1e6; 
Ls_SEP = Ls_SEP * 1e6; 
Ls_RV  = Ls_RV  * 1e6; 

V_RV = V_LV; 

%% Outputs

f(1) = abs(-V_LV - 0.5 * Vw_LV - 0.5 * Vw_SEP + Vm_SEP - Vm_LV)*1e6; 
f(2) = abs(V_RV + 0.5 * Vw_RV + 0.5 * Vw_SEP + Vm_SEP - Vm_RV)*1e6;

% For end-diastole
if flag == 0 
    f(3) = abs(Ls_LV  - SL);
    f(4) = abs(Ls_SEP - SL); 
    f(5) = abs(Ls_RV  - SL); 
end 

end 
function init = initialconditions(pars,data,eta_Vtot) 

%{ 
    This function approximates steady-state initial conditions. 
    Inputs: 
    pars        - vector of parameters 
    data        - input data structure with data and global parameters 
    Outputs: 
    init        - vector of initial conditions 
%} 


%% Fixed parameters 

fixpars = data.fixpars; 

% Minimal pressures (kPa) 
P_LV_m = fixpars(1); 
P_SA_m = fixpars(2); 
P_SV_m = fixpars(3);
P_RV_m = fixpars(4); 
P_PA_m = fixpars(5); 
P_PV_m = fixpars(6); 

% Minimal elastances (kPa m^(-3)) 
E_LV_m = fixpars(17);
E_RV_m = fixpars(18); 

%% Adjustable parameters

C_SA = pars(1); 
C_SV = pars(2); 
C_PA = pars(3); 
C_PV = pars(4); 

%% Initial conditions 

% Displacements (m)                 1 - 4 
xm_LV_0  = data.deformation.xm_LV_0;
xm_SEP_0 = data.deformation.xm_SEP_0;
xm_RV_0  = data.deformation.xm_RV_0;
ym_0     = data.deformation.ym_0; 

% Sarcomere lengths (m)             5 - 7 
Lsc_LV_0  = 2 * 1e-6;
Lsc_SEP_0 = 2 * 1e-6;
Lsc_RV_0  = 2 * 1e-6;

% Volumes (m^3)                     8 - 15 
V_LV_0 = eta_Vtot*P_LV_m / E_LV_m;
V_SA_0 = eta_Vtot*C_SA * P_SA_m;  
V_SV_0 = eta_Vtot*C_SV * P_SV_m;  
V_RV_0 = eta_Vtot*P_RV_m / E_RV_m; 
V_PA_0 = eta_Vtot*C_PA * P_PA_m;  
V_PV_0 = eta_Vtot*C_PV * P_PV_m;   

% Create initial conditions vector 
init = [xm_LV_0; xm_SEP_0; xm_RV_0; ym_0;           % 1-4
    Lsc_LV_0; Lsc_SEP_0; Lsc_RV_0;                 % 5-7
    V_LV_0; V_SA_0; V_SV_0; V_RV_0; V_PA_0; V_PV_0;   % 8-13
    ]; 

% Solve system of Triseg equations to determine consistent initialization
% for DAE 
x0   = log(init(1:7)); 
opts = optimoptions('fsolve','Display','off',...
    'MaxFunctionEvaluations',2e3); 
xopt = fsolve(@(x) triseg(x,pars,data,init),x0,opts); 
init(1:7) = exp(xopt(1:7)); 

end 
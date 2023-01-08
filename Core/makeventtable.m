function venttable = makeventtable(loading,data)

%{ 

This function creates a table with all of the important metrics for 
both ventricles calculated from the model. 

Inputs: 
outputs     - output structure from model_sol.m 
data        - input data structure with data and global parameters 
    
Outputs: 
venttable   - table with calculated metrics for both ventricles. 

%}

%% Unpack data structure

HR = data.HR; 

%% Unpack loading structure 
base_loc = loading.base_loc; 

i_ES_vec = loading.i_ES_vec; i_ES = i_ES_vec(base_loc); 
i_ED_vec = loading.i_ED_vec; i_ED = i_ED_vec(base_loc); 

V_LV_mat = loading.V_LV_mat;
V_RV_mat = loading.V_RV_mat; 

P_LV_mat = loading.P_LV_mat; 
P_RV_mat = loading.P_RV_mat; 

SV_LV_vec = loading.SV_LV_vec; 
SV_RV_vec = loading.SV_RV_vec; 

EF_LV_vec = loading.EF_LV_vec; 
EF_RV_vec = loading.EF_RV_vec; 

CO_LV_vec = loading.CO_LV_vec; 
CO_RV_vec = loading.CO_RV_vec;

CP_LV_vec = loading.CP_LV_vec; 
CP_RV_vec = loading.CP_RV_vec; 

Cm_SEP_mat = loading.Cm_SEP_mat;

EDV_LV_vec = loading.EDV_LV_vec; 
EDP_LV_vec = loading.EDP_LV_vec; 
EDV_RV_vec = loading.EDV_RV_vec; 
EDP_RV_vec = loading.EDP_RV_vec; 

ESV_LV_vec = loading.ESV_LV_vec; 
ESP_LV_vec = loading.ESP_LV_vec; 
ESV_RV_vec = loading.ESV_RV_vec; 
ESP_RV_vec = loading.ESP_RV_vec; 

% Volumes (mL)
V_LV = V_LV_mat(:,base_loc);
V_RV = V_RV_mat(:,base_loc);
    
% Pressures (mmHg)
P_LV = P_LV_mat(:,base_loc);  
P_RV = P_RV_mat(:,base_loc);
    
% Stroke volume (mL) 
SV_LV = SV_LV_vec(base_loc);
SV_RV = SV_RV_vec(base_loc);
    
% Ejection fraction (dimensionless) 
EF_LV = EF_LV_vec(base_loc);
EF_RV = EF_RV_vec(base_loc);
    
% Cardiac output (L min^(-1)) 
CO_LV = CO_LV_vec(base_loc);  
CO_RV = CO_RV_vec(base_loc);     

% Cardiac power (W)  
CP_LV = CP_LV_vec(base_loc); 
CP_RV = CP_RV_vec(base_loc); 

% Septal Curvature 
Cm_SEP_vec = Cm_SEP_mat(:,base_loc); 
Cm_SEP_min = min(Cm_SEP_vec);

% EDV, EDP, ESV, ESP
EDV_LV = EDV_LV_vec(base_loc); 
EDP_LV = EDP_LV_vec(base_loc); 
EDV_RV = EDV_RV_vec(base_loc); 
EDP_RV = EDP_RV_vec(base_loc); 

ESV_LV = ESV_LV_vec(base_loc); 
ESP_LV = ESP_LV_vec(base_loc); 
ESV_RV = ESV_RV_vec(base_loc); 
ESP_RV = ESP_RV_vec(base_loc); 

    
%% Make output table 

Vent = ['LV';'RV']; 
SV = [SV_LV; SV_RV]; 
EF = [EF_LV; EF_RV];
CO = [CO_LV; CO_RV]; 
CP = [CP_LV; CP_RV]; 
ESP = [ESP_LV; ESP_RV];
ESV = [ESV_LV; ESV_RV];
EDP = [EDP_LV; EDP_RV];
EDV = [EDV_LV; EDV_RV];

Cm_SEP = [Cm_SEP_min; 0];

% [EDV, EDP, ESV, ESP] = getEDESvals(V_LV,V_RV,P_LV,P_RV,i_ES,i_ED); 

venttable = table(Vent,SV,EF,CO,CP,ESP,EDP,ESV,EDV,Cm_SEP);

end 

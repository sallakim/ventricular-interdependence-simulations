function loading = volumeloading(pars,data,eta_Vtot_vec)

%{
    This function executes the volume loading experiment for the healthy
    and dysfunction cases. This function incorporates running the model in
    parallel. 
    Inputs: 
    pars            - parameter vector 
    data            - input data structure with data and global parameters 
    eta_Vtot_vec    - vector of scalar to progressively increase total
                      blood volume 
    Outputs: 
    loading         - output structure for volume loading experiments 
%}

%% Run volume loading experiment 

% Set time vector for up to 2 seconds 
time = 0:data.dt:2; 

HR = data.HR;

T = 60/HR; 

% Initialize vectors 
EDV_LV_vec = zeros(size(eta_Vtot_vec)); EDV_RV_vec = zeros(size(eta_Vtot_vec)); 
EDP_LV_vec = zeros(size(eta_Vtot_vec)); EDP_RV_vec = zeros(size(eta_Vtot_vec)); 
ESV_LV_vec = zeros(size(eta_Vtot_vec)); ESV_RV_vec = zeros(size(eta_Vtot_vec)); 
ESP_LV_vec = zeros(size(eta_Vtot_vec)); ESP_RV_vec = zeros(size(eta_Vtot_vec)); 
    
SV_LV_vec = zeros(size(eta_Vtot_vec)); 
SV_RV_vec = zeros(size(eta_Vtot_vec)); 

Vtot_vec = zeros(size(eta_Vtot_vec)); 
P_SAbar_vec = zeros(size(eta_Vtot_vec)); 

V_LV_mat = zeros(length(time),length(eta_Vtot_vec)); 
V_RV_mat = zeros(length(time),length(eta_Vtot_vec)); 
P_LV_mat = zeros(length(time),length(eta_Vtot_vec)); 
P_RV_mat = zeros(length(time),length(eta_Vtot_vec)); 

V_PA_mat = zeros(length(time),length(eta_Vtot_vec)); 
V_PV_mat = zeros(length(time),length(eta_Vtot_vec)); 

P_SA_mat = zeros(length(time),length(eta_Vtot_vec)); 
P_SV_mat = zeros(length(time),length(eta_Vtot_vec)); 
P_PA_mat = zeros(length(time),length(eta_Vtot_vec)); 
P_PV_mat = zeros(length(time),length(eta_Vtot_vec)); 
P_peri_mat = zeros(length(time),length(eta_Vtot_vec)); 

Q_m_valve_mat = zeros(length(time),length(eta_Vtot_vec)); 
Q_a_valve_mat = zeros(length(time),length(eta_Vtot_vec)); 
Q_t_valve_mat = zeros(length(time),length(eta_Vtot_vec)); 
Q_p_valve_mat = zeros(length(time),length(eta_Vtot_vec)); 

i_ES_vec = zeros(size(eta_Vtot_vec));
i_ED_vec = zeros(size(eta_Vtot_vec));

CO_LV_vec = zeros(size(eta_Vtot_vec));
CO_RV_vec = zeros(size(eta_Vtot_vec));

CP_LV_vec = zeros(size(eta_Vtot_vec));
CP_RV_vec = zeros(size(eta_Vtot_vec));

Cm_SEP_mat = zeros(length(time),length(eta_Vtot_vec)); 
Cm_LV_mat = zeros(length(time),length(eta_Vtot_vec)); 
Cm_RV_mat = zeros(length(time),length(eta_Vtot_vec)); 


for i = 1:length(eta_Vtot_vec) 

    % Get parameters 
    eta_Vtot = eta_Vtot_vec(i); 

    % Solve model 
    outputs_loading = volumeloading_wrap(pars,data,eta_Vtot);

    time = outputs_loading.time;

    V_LV = outputs_loading.volumes.V_LV;    V_RV = outputs_loading.volumes.V_RV; 
    P_LV = outputs_loading.pressures.P_LV;  P_RV = outputs_loading.pressures.P_RV; 

    V_LV_mat(:,i) = V_LV;       V_RV_mat(:,i) = V_RV; 
    P_LV_mat(:,i) = P_LV;       P_RV_mat(:,i) = P_RV; 
    
    V_PA_mat(:,i) = outputs_loading.volumes.V_PA; 
    V_PV_mat(:,i) = outputs_loading.volumes.V_PV; 

    P_SA = outputs_loading.pressures.P_SA;    P_PA = outputs_loading.pressures.P_PA;
    P_SV = outputs_loading.pressures.P_SV;    P_PV = outputs_loading.pressures.P_PV;
    P_peri = outputs_loading.pressures.P_peri; 

    P_SA_mat(:,i) = P_SA;       P_PA_mat(:,i) = P_PA; 
    P_SV_mat(:,i) = P_SV;       P_PV_mat(:,i) = P_PV; 
    P_peri_mat(:,i) = P_peri; 

    Q_m_valve = outputs_loading.flows.Q_m_valve; 
    Q_a_valve = outputs_loading.flows.Q_a_valve; 
    Q_t_valve = outputs_loading.flows.Q_t_valve; 
    Q_p_valve = outputs_loading.flows.Q_p_valve; 

    Q_m_valve_mat(:,i) = Q_m_valve; 
    Q_a_valve_mat(:,i) = Q_a_valve; 
    Q_t_valve_mat(:,i) = Q_t_valve; 
    Q_p_valve_mat(:,i) = Q_p_valve; 

%     i_ES = find(diff(Q_m_valve) > 0,1,'first'); 
%     i_ED = find(diff(Q_a_valve) > 0,1,'first'); 
%     Q_a_valve_x = Q_a_valve(i_ED+1000:end);
%     x_vec = ones(1,i_ED + 1000-1);
%     Q_a_valve_x = [x_vec,Q_a_valve_x];
%     i_ED = find(diff(Q_a_valve_x) > 0,1,'first');  
     
%     ED = time(i_ED); 
%     ED_2 = ED + T; 
%     i_ED = find(time == ED_2);

    i_ES_LV = find(diff(Q_a_valve) < 0,1,'last');
    i_ED_LV = find(diff(Q_m_valve) < 0,1,'last');

    i_ES_RV = find(diff(Q_p_valve) < 0,1,'last'); %last
    i_ED_RV = find(diff(Q_t_valve) < 0,1,'last'); %first
%     
    EDV_LV = V_LV(i_ED_LV);        EDV_RV = V_RV(i_ED_RV); 
    EDP_LV = P_LV(i_ED_LV);        EDP_RV = P_RV(i_ED_RV);
    ESV_LV = V_LV(i_ES_LV);        ESV_RV = V_RV(i_ES_RV);
    ESP_LV = P_LV(i_ES_LV);        ESP_RV = P_RV(i_ES_RV);
%     
%      [EDV, EDP, ESV, ESP] = getEDESvals(V_LV,V_RV,P_LV,P_RV);
% 
%     EDV_LV = EDV(1,:);          EDV_RV = EDV(2,:); 
%     EDP_LV = EDP(1,:);          EDP_RV = EDP(2,:); 
%     ESV_LV = ESV(1,:);          ESV_RV = ESV(2,:); 
%     ESP_LV = ESP(1,:);          ESP_RV = ESP(2,:);
    
    EDV_LV_vec(i) = EDV_LV;   EDV_RV_vec(i) = EDV_RV; 
    EDP_LV_vec(i) = EDP_LV;   EDP_RV_vec(i) = EDP_RV; 
    ESV_LV_vec(i) = ESV_LV;   ESV_RV_vec(i) = ESV_RV; 
    ESP_LV_vec(i) = ESP_LV;   ESP_RV_vec(i) = ESP_RV; 

    SV_LV_vec(i) = max(V_LV) - min(V_LV); 
    SV_RV_vec(i) = max(V_RV) - min(V_RV); 

    Vtot_vec(i,:) = outputs_loading.volumes.Vtot; 

    P_SA_M = max(outputs_loading.pressures.P_SA); 
    P_SA_m = min(outputs_loading.pressures.P_SA); 
    P_SAbar_vec(i) = (1/3) * P_SA_M + (2/3) * P_SA_m; % To compare with Pbar

    Q_a_valve = outputs_loading.flows.Q_a_valve;    
    Q_p_valve = outputs_loading.flows.Q_p_valve; 

    CO_LV_vec(i) = trapz(time/60,Q_a_valve)/(time(end)/60 - time(1)/60); %SV * HR_end * 1e-3 
    CO_RV_vec(i) = trapz(time/60,Q_p_valve)/(time(end)/60 - time(1)/60); %SV * HR_end * 1e-3 
    
    CP_LV = trapz(P_LV,V_LV) / 7.5 * 1e-3 * HR/60; %mean(P_sa(beat)) / 7.5 * 1e3 * SV * 1e-6 * HR_end/60 
    CP_RV = trapz(P_RV,V_RV) / 7.5 * 1e-3 * HR/60;

    CP_LV = CP_LV/2; % average over two beats
    CP_RV = CP_RV/2; 
    
    CP_LV_vec(i) = CP_LV; 
    CP_RV_vec(i) = CP_RV; 

    Cm_SEP = outputs_loading.curvatures.Cm_SEP;
    Cm_LV  = outputs_loading.curvatures.Cm_LV; 
    Cm_RV  = outputs_loading.curvatures.Cm_RV; 

    Cm_SEP_mat(:,i) = Cm_SEP;
    Cm_LV_mat(:,i) = Cm_LV;
    Cm_RV_mat(:,i) = Cm_RV;

end

y_v = outputs_loading.activation.y_v;

%% Find run that has a mean pressure closest to healthy case 

SAbar = 95; 
[~,m] = min(abs(P_SAbar_vec - SAbar)); 

EF_LV = SV_LV_vec./EDV_LV_vec; 
EF_RV = SV_RV_vec./EDV_RV_vec; 


%% Create structure from volume loading experiments 

loading.base_loc = m; 

loading.EDV_LV_vec = EDV_LV_vec; 
loading.EDV_RV_vec = EDV_RV_vec; 
loading.EDP_LV_vec = EDP_LV_vec; 
loading.EDP_RV_vec = EDP_RV_vec; 
loading.ESV_LV_vec = ESV_LV_vec; 
loading.ESV_RV_vec = ESV_RV_vec; 
loading.ESP_LV_vec = ESP_LV_vec; 
loading.ESP_RV_vec = ESP_RV_vec; 
    
loading.SV_LV_vec = SV_LV_vec;
loading.SV_RV_vec = SV_RV_vec; 
loading.Vtot_vec = Vtot_vec; 
    
loading.V_LV_mat = V_LV_mat; 
loading.V_RV_mat = V_RV_mat; 
loading.P_LV_mat = P_LV_mat; 
loading.P_RV_mat = P_RV_mat; 

loading.V_PA_mat = V_PA_mat;
loading.V_PV_mat = V_PV_mat;

loading.P_SA_mat = P_SA_mat; 
loading.P_SV_mat = P_SV_mat; 
loading.P_PA_mat = P_PA_mat; 
loading.P_PV_mat = P_PV_mat; 
loading.P_peri_mat = P_peri_mat; 

loading.Q_m_valve_mat = Q_m_valve_mat;
loading.Q_a_valve_mat = Q_a_valve_mat;
loading.Q_t_valve_mat = Q_t_valve_mat; 
loading.Q_p_valve_mat = Q_p_valve_mat; 

loading.EF_LV_vec = EF_LV;
loading.EF_RV_vec = EF_RV; 

loading.CO_LV_vec = CO_LV_vec; 
loading.CO_RV_vec = CO_RV_vec; 

loading.CP_LV_vec = CP_LV_vec;
loading.CP_RV_vec = CP_RV_vec;

loading.Cm_SEP_mat = Cm_SEP_mat;
loading.Cm_LV_mat  = Cm_LV_mat;
loading.Cm_RV_mat  = Cm_RV_mat;

loading.time = time;

loading.T = data.T; 

loading.y_v = y_v; 

loading.i_ES_vec = i_ES_vec; 
loading.i_ED_vec = i_ED_vec;

% output for visualization
loading.displacements = outputs_loading.displacements;
loading.tensions = outputs_loading.tensions;
end 



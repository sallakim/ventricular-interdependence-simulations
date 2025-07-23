%{ 
RUN_invivo_VL

Driver to run the volume loading experiment for the in vivo version of the 
cardiac and cardiovascular model. This code requires the user to choose
which experiment to run, which can be done in the next section. 
%}

clear all
% close all

tic 

addpath Core
addpath InVivo

n_cores = 4; 

%% Choose experiment to run

% run_experiment = 'Healthy';
% run_experiment = 'LVSD'; 
% run_experiment = 'LVDD'; 
% run_experiment = 'RVSD'; 
% run_experiment = 'RVDD'; 
run_experiments = {'Healthy', 'LVSD', 'LVDD', 'RVSD', 'RVDD'}; 
% run_experiments = {'Healthy'}; 


delete(gcp('nocreate'))
% p = parpool('local',n_cores);

%% Flags 

printoutfigs_on = 1; % 0 - off; 1 - on 

for i_exp = 1:length(run_experiments)
%% Repeat for all
run_experiment = run_experiments{i_exp};

%% Load Data 

data = makedatastructure;

data.printoutfigs_on = printoutfigs_on; 

%% Global parameters

ODE_TOL = 1e-8; 
data.gpars.ODE_TOL = ODE_TOL; 
eta_Vtot_vec = [1:0.1:3.5];
%% Volume loading for disease case 
if ~strcmp('Healthy',run_experiment)
    % Redistribute blood volume for disease case 
    bvd_SA = .6; 
    bvd_PA = .35;
    bvd_SV = .8;
    bvd_PV = .75;

    data.bvd_SA = bvd_SA; 
    data.bvd_PA                                                                     = bvd_PA; 
    data.bvd_SV = bvd_SV; 
    data.bvd_PV = bvd_PV; 

    data_m = data; 
    data_s = data;
end

switch run_experiment
    case 'Healthy'
        %% Volume loading for healthy case
        % eta_Vtot_vec = 1; 

        [pars_h,~,~,data_h] = parameters(data);
        loading_h = volumeloading(pars_h,data_h,eta_Vtot_vec);

        TriSegAnimExport(loading_h, 'normal');

    case 'LVSD'
        % Moderate
        eta_k_act_LV = 0.6;
        data_m.eta_k_act_LV = eta_k_act_LV;
        [pars_m,~,~,data_m] = parameters(data_m); 
        loading_m = volumeloading(pars_m,data_m,eta_Vtot_vec);
        
        % Severe 
        eta_k_act_LV = 0.4;
        data_s.eta_k_act_LV = eta_k_act_LV;
        [pars_s,~,~,data_s] = parameters(data_s); 
        loading_s = volumeloading(pars_s,data_s,eta_Vtot_vec);
        TriSegAnimExport(loading_s, 'LVSD');

    case 'LVDD'
        % Moderate
        eta_k_pas_LV = 7;
        data_m.eta_k_pas_LV = eta_k_pas_LV;
        [pars_m,~,~,data_m] = parameters(data_m); 
        loading_m = volumeloading(pars_m,data_m,eta_Vtot_vec);
        
        % Severe
        eta_k_pas_LV = 15;
        data_s.eta_k_pas_LV = eta_k_pas_LV;
        [pars_s,~,~,data_s] = parameters(data_s); 
        loading_s = volumeloading(pars_s,data_s,eta_Vtot_vec);
        TriSegAnimExport(loading_s, 'LVDD');

    case 'RVSD'
        
        % Moderate
        eta_k_act_RV = 0.6;
        data_m.eta_k_act_RV = eta_k_act_RV;
        [pars_m,~,~,data_m] = parameters(data_m); 
        loading_m = volumeloading(pars_m,data_m,eta_Vtot_vec);
        
        % Severe
        eta_k_act_RV = 0.4 ;
        data_s.eta_k_act_RV = eta_k_act_RV;
        [pars_s,~,~,data_s] = parameters(data_s); 
        loading_s = volumeloading(pars_s,data_s,eta_Vtot_vec); 
        TriSegAnimExport(loading_s, 'RVSD');

    case 'RVDD'
        % Moderate
        eta_k_pas_RV = 7;
        data_m.eta_k_pas_RV = eta_k_pas_RV;
        [pars_m,~,~,data_m] = parameters(data_m); 
        loading_m = volumeloading(pars_m,data_m,eta_Vtot_vec);
        
        % Severe
        eta_k_pas_RV = 15;
        data_s.eta_k_pas_RV = eta_k_pas_RV;
        [pars_s,~,~,data_s] = parameters(data_s); 
        loading_s = volumeloading(pars_s,data_s,eta_Vtot_vec);
        TriSegAnimExport(loading_s, 'RVDD');

end

delete(gcp('nocreate'))


%% Plot volume loading results 

if strcmp('Healthy',run_experiment)
    venttable_h = makeventtable(loading_h,data_h)

    loading = loading_h; 
    loading.run_experiment = run_experiment; 
    loading.data = data_h; 

    plot_invivo_VL_healthyonly(loading)

else 
    venttable_h = makeventtable(loading_h,data_h)
    venttable_m = makeventtable(loading_m,data_m)
    venttable_s = makeventtable(loading_s,data_s)

    loading.run_experiment = run_experiment; 
    loading.healthy  = loading_h; 
    loading.moderate = loading_m; 
    loading.severe   = loading_s; 
    
    loading.pars_h = pars_h; 
    loading.pars_m = pars_m; 
    loading.pars_s = pars_s; 
    
    data.data_h = data_h; 
    data.data_m = data_m; 
    data.data_s = data_s; 
    
    plot_invivo_VL_healthymodsev(loading,data)
end

elapsed_time = toc/60
end
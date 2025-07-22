% Modified TriSeg Model 
% driver: This code is the main driver for the cardiovascular mechanics model 

clear all
%close all

addpath ExVivo
addpath Core

printon = 1; 

%% Load pseudodata

data = makedatastructure; 

EDV_LV = data.EDV_LV; 
EDV_RV = data.EDV_RV; 
EDP_LV = data.EDP_LV; 
EDP_RV = data.EDP_RV;

%% Global parameters

ODE_TOL = 1e-8; 
data.gpars.ODE_TOL = ODE_TOL; 


%% Solve model 

V_LV_vec = []; 
V_RV_vec = []; 
P_LV_vec = []; 
P_RV_vec = []; 

eta_EDV = [.5 1 1.5 2];
eta_Vtot = [.75:.125:1.5]; 
for q = 1:length(eta_EDV)
    
    data.EDV_LV = eta_EDV(q)*EDV_LV;
    data.EDV_RV = eta_EDV(q)*EDV_RV;
    
    %% Get Parameters 

    [adjpars,~,~,data] = parameters(data); 

for i = 1:length(eta_Vtot)
    data.eta_Vtot = eta_Vtot(i); 

    outputs = model_sol_exvivo(adjpars,data); 
    
    VLV1 = outputs.volumes.V_LV;
    PLV1 = outputs.pressures.P_LV; 

    VRV1 = outputs.volumes.V_RV;
    PRV1 = outputs.pressures.P_RV; 

    V_LV_vec(i) = VLV1; 
    V_RV_vec(i) = VRV1; 
    P_LV_vec(i) = PLV1; 
    P_RV_vec(i) = PRV1; 

end

norm_loc = find(eta_Vtot == 1);
An = 28; 
Bn = 3; 

% LV
EDP_LV_norm = P_LV_vec(norm_loc); 
EDV_LV_norm = V_LV_vec(norm_loc);

V_0  = EDV_LV_norm * (0.6 - 0.006 * EDP_LV_norm); 
V_30 = V_0 + (EDV_LV_norm - V_0) / ((EDP_LV_norm / An)^(1/Bn));
a_EDV_LV_normalized = (V_LV_vec - V_0)./(V_30 - V_0);

% RV
EDP_RV_norm = P_RV_vec(norm_loc); 
EDV_RV_norm = V_RV_vec(norm_loc);

V_0  = EDV_RV_norm * (0.6 - 0.006 * EDP_RV_norm); 
V_30 = V_0 + (EDV_RV_norm - V_0) / ((EDP_RV_norm / An)^(1/Bn)); 
a_EDV_RV_normalized = (V_RV_vec - V_0)./(V_30 - V_0);

save (['EDV' num2str(100*eta_EDV(q))])

end

%% Klotz curve

An = 28; 
Bn = 3; 

% LV
EDP_LV_norm = EDP_LV; 
EDV_LV_norm = EDV_LV; 

V_0  = EDV_LV_norm * (0.6 - 0.006 * EDP_LV_norm); 
V_30 = V_0 + (EDV_LV_norm - V_0) / ((EDP_LV_norm / An)^(1/Bn));
Beta_lv = log(EDP_LV_norm/30) / log(EDV_LV_norm / V_30); 
Alpha_lv = 30 / V_30^Beta_lv; 

% RV
EDP_RV_normal = EDP_RV; 
EDV_RV_normal = EDV_RV; 

V_0  = EDV_RV_normal * (0.6 - 0.006 * EDP_RV_normal); 
V_30 = V_0 + (EDV_RV_normal - V_0) / ((EDP_RV_normal / An)^(1/Bn));
Beta_RV = log(EDP_RV_normal/30) / log(EDV_RV_normal / V_30); 
Alpha_RV = 30 / V_30^Beta_RV; 

V_EDPVR = [50:230]; 
P_LV_EDPVR = Alpha_lv * V_EDPVR.^Beta_lv;
P_RV_EDPVR = Alpha_RV * V_EDPVR.^Beta_RV;

%% 
hfig100 = figure(100);
clf
hold on 
load EDV50.mat
plot(V_LV_vec,P_LV_vec,'r--')
load EDV100.mat
plot(V_LV_vec,P_LV_vec,'r-x')
load EDV150.mat
plot(V_LV_vec,P_LV_vec,'r-*')
load EDV200.mat
plot(V_LV_vec,P_LV_vec,'r-o')
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
ylim([0 50])
xlim([0 350])
legend('x0.5','x1.0','x1.5','x2.0','location','northwest')
set(gca,'FontSize',20)

hfig101 = figure(101);
clf
hold on 
load EDV50.mat
plot(V_RV_vec,P_RV_vec,'b--')
load EDV100.mat
plot(V_RV_vec,P_RV_vec,'b-x')
load EDV150.mat
plot(V_LV_vec,P_RV_vec,'b-*')
load EDV200.mat
plot(V_LV_vec,P_RV_vec,'b-o')
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
ylim([0 50])
xlim([0 350])
legend('x0.5','x1.0','x1.5','x2.0','location','northwest')
set(gca,'FontSize',20)

figure102 = figure(102);
clf
hold on 
load EDV50.mat
plot(a_EDV_LV_normalized,P_LV_vec,'r--')
load EDV100.mat
plot(a_EDV_LV_normalized,P_LV_vec,'r-x')
load EDV150.mat
plot(a_EDV_LV_normalized,P_LV_vec,'r-*')
load EDV200.mat
plot(a_EDV_LV_normalized,P_LV_vec,'r-o')
xlabel('Normalized Volume (mL)')
ylabel('Pressure (mmHg)')
ylim([0 50])
xlim([0 1.2])
legend('x0.5','x1.0','x1.5','x2.0','location','northwest')
set(gca,'FontSize',20)

hfig103 = figure(103);
clf
hold on 
load EDV50.mat
plot(a_EDV_RV_normalized,P_RV_vec,'b--')
load EDV100.mat
plot(a_EDV_RV_normalized,P_RV_vec,'b-x')
load EDV150.mat
plot(a_EDV_RV_normalized,P_RV_vec,'b-*')
load EDV200.mat
plot(a_EDV_RV_normalized,P_RV_vec,'b-o')
xlabel('Normalized Volume (mL)')
ylabel('Pressure (mmHg)')
ylim([0 50])
xlim([0 1.2])
legend('x0.5','x1.0','x1.5','x2.0','location','northwest')
set(gca,'FontSize',20)

 print(hfig100,'-dpng',strcat('Figures/','ExVivo','/Figure2_c.png'))
 print(hfig101,'-dpng',strcat('Figures/','ExVivo','/Figure2_d.png'))
 print(hfig102,'-dpng',strcat('Figures/','ExVivo','/Figure2_e.png'))
 print(hfig103,'-dpng',strcat('Figures/','ExVivo','/Figure2_f.png'))


 print(hfig100,'-depsc2',strcat('Figures/','ExVivo','/Figure2_c.eps'))
 print(hfig101,'-depsc2',strcat('Figures/','ExVivo','/Figure2_d.eps'))
 print(hfig102,'-depsc2',strcat('Figures/','ExVivo','/Figure2_e.eps'))
 print(hfig103,'-depsc2',strcat('Figures/','ExVivo','/Figure2_f.eps'))

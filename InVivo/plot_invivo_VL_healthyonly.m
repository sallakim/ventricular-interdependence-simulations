function [] = plot_invivo_VL_healthyonly(loading)

%{ 
    Plot the volume loading experiment for the in vivo version of the model
    for the healthy case only.
    Inputs: 
    loading         - input structure from volumeloading.m containing 
                      vectors for clinical metrics 
    data            - input data structure with data and global parameters 
    Outputs: (figures) 
%} 

%% Flags

printoutfigs_on = loading.data.printoutfigs_on; 

%% Unpack data

SPbar = loading.data.SPbar*7.5;
DPbar = loading.data.DPbar*7.5; 

T = 60/loading.data.HR; 

%% Unpack healthy loading structure 

run_experiment = loading.run_experiment; 

time = loading.time;

base_loc = loading.base_loc; 

EDV_LV_vec = loading.EDV_LV_vec;    EDV_RV_vec = loading.EDV_RV_vec; 
EDP_LV_vec = loading.EDP_LV_vec;    EDP_RV_vec = loading.EDP_RV_vec; 
ESV_LV_vec = loading.ESV_LV_vec;    ESV_RV_vec = loading.ESV_RV_vec; 
ESP_LV_vec = loading.ESP_LV_vec;    ESP_RV_vec = loading.ESP_RV_vec; 

SV_LV_vec = loading.SV_LV_vec;      SV_RV_vec = loading.SV_RV_vec; 

V_LV_mat = loading.V_LV_mat;        V_RV_mat = loading.V_RV_mat;
P_LV_mat = loading.P_LV_mat;        P_RV_mat = loading.P_RV_mat;
Vtot_vec = loading.Vtot_vec;

V_PA_mat = loading.V_PA_mat;        V_PV_mat = loading.V_PV_mat;
V_PA = V_PA_mat(:,base_loc);        V_PV = V_PV_mat(:,base_loc);

P_SA_mat = loading.P_SA_mat;        P_PA_mat = loading.P_PA_mat; 
P_SV_mat = loading.P_SV_mat;        P_PV_mat = loading.P_PV_mat; 

Q_m_valve_mat = loading.Q_m_valve_mat;
Q_a_valve_mat = loading.Q_a_valve_mat;
Q_t_valve_mat = loading.Q_t_valve_mat;
Q_p_valve_mat = loading.Q_p_valve_mat;

V_LV = V_LV_mat(:,base_loc);  V_RV = V_RV_mat(:,base_loc);
P_LV = P_LV_mat(:,base_loc);  P_RV = P_RV_mat(:,base_loc);

P_SA = P_SA_mat(:,base_loc);  P_PA = P_PA_mat(:,base_loc); 
P_SV = P_SV_mat(:,base_loc);  P_PV = P_PV_mat(:,base_loc); 

Q_m_valve = Q_m_valve_mat(:,base_loc);
Q_a_valve = Q_a_valve_mat(:,base_loc);
Q_t_valve = Q_t_valve_mat(:,base_loc);
Q_p_valve = Q_p_valve_mat(:,base_loc);

CO_LV_vec = loading.CO_LV_vec;  CO_RV_vec = loading.CO_RV_vec; 
CO_LV = CO_LV_vec(base_loc);    CO_RV = CO_RV_vec(base_loc); 

CP_LV_vec = loading.CP_LV_vec;  CP_RV_vec = loading.CP_RV_vec; 
CP_LV = CP_LV_vec(base_loc);    CP_RV = CP_RV_vec(base_loc);

EF_LV_vec = loading.EF_LV_vec;  EF_RV_vec = loading.EF_RV_vec;
EF_LV = EF_LV_vec(base_loc);    EF_RV = EF_RV_vec(base_loc);

Cm_SEP_mat = loading.Cm_SEP_mat; 

Cm_SEP_vec = Cm_SEP_mat(:,base_loc);

y_v = loading.y_v; 

% Find ejection
i_AVO = find(diff(Q_a_valve) > 0,1,'first'); 
i_AVC = find(diff(Q_a_valve) < 0,1,'last');

T_AVO = time(i_AVO); 
T_AVC = time(i_AVC) - T; 

%% Make EDPVR Klotz curves

V_EDPVR = [10:200];

% [EDV_h,EDP_h] = getEDESvals(V_LV,V_RV,P_LV,P_RV,i_ES,i_ED); 

EDV_LV_h = EDV_LV_vec(base_loc);    EDV_RV_h = EDV_RV_vec(base_loc);
EDP_LV_h = EDP_LV_vec(base_loc);    EDP_RV_h = EDP_RV_vec(base_loc);
EDV_h = [EDV_LV_h, EDV_RV_h];
EDP_h = [EDP_LV_h, EDP_RV_h];

[P_LV_EDPVR_h,P_RV_EDPVR_h] = makeKlotzcurve(EDV_h,EDP_h,V_EDPVR);

%% ESPVR curves

xfit = linspace(30,100); 
    
p_LV_h = polyfit(ESV_LV_vec(2:end),ESP_LV_vec(2:end),1);
p_RV_h = polyfit(ESV_RV_vec(2:end),ESP_RV_vec(2:end),1);
yfit_LV_h = p_LV_h(1)*xfit + p_LV_h(2); 
yfit_RV_h = p_RV_h(1)*xfit + p_RV_h(2);     

%% Plot figures

gray = [.5 .5 .5]; 
green = [0 .5 0]; 
purple= [0.4940 0.1840 0.5560];

%% PV loops
hfig1 = figure(1); 
clf
hold on 
plot(V_EDPVR,P_LV_EDPVR_h,'color',gray,'linewidth',1)
plot(V_EDPVR,P_RV_EDPVR_h,'color',gray,'linewidth',1)
h1 = plot(V_LV, P_LV, 'r',   'linewidth',2);
h2 = plot(V_RV, P_RV, 'b',   'linewidth',2);
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
legend([h1 h2],'LV','RV','location','Northwest')
xlim([0 200])
ylim([0 150])
set(gca,'FontSize',20)

%% Volume traces 
hfig2 = figure(2); 
clf
hold on 
h1 = plot(loading.time,V_LV,'r-','linewidth',2);
h2 = plot(loading.time,V_RV,'b-','linewidth',2);


legend([h1,h2],'LV','RV')

xlabel('Time (s)')
ylabel('Volume (mL)')
set(gca,'FontSize',20)


%% High pressures
% LV
hfig3 = figure(3); 
clf
hold on 
h1 = plot(time,P_LV,'r-','linewidth',2);
h2 = plot(time,P_SA,'LineWidth',2,'color',"#EDB120");

legend([h1,h2],'P_{LV}','P_{SA}')

xlabel('Time (s)')
ylabel('Pressure (mmHg)')

set(gca,'FontSize',20)
ylim([0 125])


%% Low pressures 
hfig4 = figure(4); 
clf
hold on 
h1 = plot(time,P_RV,'b-','linewidth',2);
h2 = plot(time,P_SV,'m-','LineWidth',2); 
h3 = plot(time,P_PA,'c-','LineWidth',2);
h4 = plot(loading.time,P_PV,'k-','LineWidth',2);


legend([h1,h2,h3,h4],'P_{RV}','P_{SV}','P_{PA}','P_{PV}')
xlabel('Time (s)')
ylabel('Pressure (mmHg)')

set(gca,'FontSize',20)

%% Flows
hfig5 = figure(5); 
clf
hold on
plot(time,Q_m_valve,'r','linewidth',2)
plot(time,Q_a_valve,'r--','linewidth',2)
plot(time,Q_t_valve,'b','linewidth',2)
plot(time,Q_p_valve,'b--','linewidth',2)

ylabel('Flow (mL s^{-1})')
xlabel('Time (s)')
legend('Q_{m}','Q_{a}','Q_{t}','Q_{p}')

set(gca,'FontSize',20)


%% Activation function 
hfig6 = figure(6); 
clf
plot(time,y_v,'g','linewidth',3)
set(gca,'FontSize',20)
xlabel('Time (s)')
ylabel('y_{v}(t)')
title('Cardiac Activation')
xlim([0 1])

%% Volume change 
hfig7 = figure(7); 
clf
hold on 
plot([1:length(Vtot_vec)], Vtot_vec,'k')
plot([1:length(Vtot_vec)], Vtot_vec,'.m','MarkerSize',10)
h_norm = plot(base_loc,Vtot_vec(base_loc),'mo','MarkerSize',10,'Linewidth',2);
legend(h_norm,'Baseline','Location','northwest')
ylabel('Total Blood Volume (L)')
xlabel('Iteration')
set(gca,'FontSize',20)

%% PV loop progression 
% LV
hfig8 = figure(8); 
clf 
hold on 
plot(V_LV_mat,P_LV_mat,'color',gray)
h_norm = plot(V_LV_mat(:,base_loc),P_LV_mat(:,base_loc),'r','Linewidth',2);
legend(h_norm,'Baseline')
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
set(gca,'FontSize',20)
xlim([min(min(V_LV_mat))-10 max(max(V_LV_mat))+10])
ylim([0 max(max(P_LV_mat))+10])

% RV

hfig9 = figure(9); 
clf 
hold on 
plot(V_RV_mat,P_RV_mat,'color',gray)
h_norm = plot(V_RV_mat(:,base_loc),P_RV_mat(:,base_loc),'b','Linewidth',2);
legend(h_norm,'Baseline')
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
set(gca,'FontSize',20)
xlim([min(min(V_RV_mat))-10 max(max(V_RV_mat))+10])
ylim([0 max(max(P_RV_mat))+10])

%% EDPVR Plots
hfig10 = figure(10);
clf
hold on 
plot(V_EDPVR,P_LV_EDPVR_h,'color',gray,'linewidth',2)
plot(V_EDPVR,P_RV_EDPVR_h,'color',gray,'linewidth',2)
h1 = plot(EDV_LV_vec,EDP_LV_vec,'ro','Markersize',6,'MarkerFaceColor','r');
h2 = plot(EDV_RV_vec,EDP_RV_vec,'bo','Markersize',6,'MarkerFaceColor','b');
l_norm_h = plot(EDV_LV_vec(base_loc),EDP_LV_vec(base_loc),'ko','Markersize',10,'linewidth',2);
r_norm_h = plot(EDV_RV_vec(base_loc),EDP_RV_vec(base_loc),'ko','Markersize',10,'linewidth',2);
legend([h1 h2 l_norm_h],'LV','RV','Baseline','Location','northwest')
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
xlim([50 250])
ylim([0 50])
set(gca,'FontSize',20)


%% ESPVR Plots
hfig11 = figure(11);
clf
hold on 
h1 = plot(ESV_LV_vec,ESP_LV_vec,'ro','Markersize',6,'MarkerFaceColor','r');
h2 = plot(ESV_RV_vec,ESP_RV_vec,'bo','Markersize',6,'MarkerFaceColor','b');
    
l_norm_h = plot(ESV_LV_vec(base_loc),ESP_LV_vec(base_loc),'ko','Markersize',10,'linewidth',2);
r_norm_h = plot(ESV_RV_vec(base_loc),ESP_RV_vec(base_loc),'ko','Markersize',10,'linewidth',2);

legend([h1 h2 l_norm_h],'LV','RV','Baseline','Location','northeast')

xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')    
    
set(gca,'FontSize',20)

%% Frank Starling
% LV 
hfig12 = figure(12);
clf
hold on
plot(EDP_LV_vec, SV_LV_vec,'o','color',green','MarkerSize',6,'MarkerFaceColor',green)
h_norm_h = plot(EDP_LV_vec(base_loc),SV_LV_vec(base_loc),'ko','Markersize',10,'linewidth',2);
legend([h_norm_h],'Baseline','Location','northwest') 
xlabel('EDP (mmHg)')
ylabel('SV (mL)')
set(gca,'FontSize',20)
xlim([0 50])

% RV
hfig13 = figure(13);
clf
hold on
plot(EDP_RV_vec, SV_RV_vec,'o','color',green','MarkerSize',6,'MarkerFaceColor',green) 
h_norm_h = plot(EDP_RV_vec(base_loc),SV_RV_vec(base_loc),'ko','Markersize',10,'linewidth',2);
legend([h_norm_h],'Baseline','Location','northwest') 
xlabel('EDP (mmHg)')
ylabel('SV (mL)')
set(gca,'FontSize',20)
xlim([0 20])

%% Septal Curvature 

hfig14 = figure(14);
clf
hold on 
plot(time,Cm_SEP_vec,'color', green,'linewidth',2)
x = [T_AVO T_AVC T_AVC T_AVO];
y1 = [0 0 .5 .5];
patch(x,y1,gray,'linestyle','none')
alpha(.05)
text(T_AVO+.02,.04,'Ejection','fontsize',20)
xlabel('Time (s)')
ylabel('Curvature (cm^{-1})')
title(strcat('Septal Curvature'))
ylim([.15 .4])
xlim([0 1])
set(gca,'FontSize',20)


%% Figure EDPVR only LV/RV

hfig15 = figure(15);
clf
hold on 
h1 = plot(V_EDPVR,P_LV_EDPVR_h,'color',gray,'linewidth',2);
h2 = plot(EDV_LV_vec,EDP_LV_vec,'ro','Markersize',6,'MarkerFaceColor','r');    
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
xlim([50 200])
ylim([0 50])
legend([h1 h2],'Klotz','Model','location','northwest')
set(gca,'FontSize',20)

hfig16 = figure(16);
clf
hold on 
h1=plot(V_EDPVR,P_RV_EDPVR_h,'color',gray,'linewidth',2);
h2 = plot(EDV_RV_vec,EDP_RV_vec,'bo','Markersize',6,'MarkerFaceColor','b');
legend([h1 h2],'Klotz','Model','location','northwest')
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
xlim([50 200])
ylim([0 50])
set(gca,'FontSize',20)


%% FIGURE 3 SUBPLOT 

hfig17 = figure(17);
clf 

% PV loops 
subplot(2,2,1)
hold on
plot(V_EDPVR,P_LV_EDPVR_h,'color',gray,'linewidth',1)
plot(V_EDPVR,P_RV_EDPVR_h,'color',gray,'linewidth',1)
h1 = plot(V_LV, P_LV, 'r',   'linewidth',2);
h2 = plot(V_RV, P_RV, 'b',   'linewidth',2);
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
legend([h1 h2],'LV','RV','location','Northeast')
xlim([0 200])
ylim([0 150])
set(gca,'FontSize',8)

% Volumes
subplot(2,2,2)
hold on 
h1 = plot(loading.time,V_LV,'r-','linewidth',2);
h2 = plot(loading.time,V_RV,'b-','linewidth',2);
legend([h1,h2],'LV','RV')
xlabel('Time (s)')
ylabel('Volume (mL)')
set(gca,'FontSize',8)

% High pressures 
subplot(2,2,3)
hold on
h1 = plot(time,P_LV,'r-','linewidth',2);
h2 = plot(time,P_SA,'LineWidth',2,'color',"#EDB120");
legend([h1,h2],'P_{LV}','P_{SA}')
xlabel('Time (s)')
ylabel('Pressure (mmHg)')
set(gca,'FontSize',8)
ylim([0 125])

% Low pressures 
subplot(2,2,4)
hold on
h1 = plot(time,P_RV,'b-','linewidth',2);
h2 = plot(time,P_SV,'m-','LineWidth',2); 
h3 = plot(time,P_PA,'c-','LineWidth',2);
h4 = plot(loading.time,P_PV,'k-','LineWidth',2);
legend([h1,h2,h3,h4],'P_{RV}','P_{SV}','P_{PA}','P_{PV}')
xlabel('Time (s)')
ylabel('Pressure (mmHg)')
set(gca,'FontSize',8)

%% FIGURE 4 SUBPLOT 

hfig18 = figure(18); 
clf

% PV loop progression, LV 
subplot(3,2,1)
hold on 
plot(V_LV_mat,P_LV_mat,'color',gray)
h_norm = plot(V_LV_mat(:,base_loc),P_LV_mat(:,base_loc),'r','Linewidth',2);
legend(h_norm,'Baseline')
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
xlim([min(min(V_LV_mat))-10 max(max(V_LV_mat))+10])
ylim([0 max(max(P_LV_mat))+10])
set(gca,'FontSize',8)

% PV loop progression, RV 
subplot(3,2,2)
hold on 
plot(V_RV_mat,P_RV_mat,'color',gray)
h_norm = plot(V_RV_mat(:,base_loc),P_RV_mat(:,base_loc),'b','Linewidth',2);
legend(h_norm,'Baseline')
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
xlim([min(min(V_RV_mat))-10 max(max(V_RV_mat))+10])
ylim([0 max(max(P_RV_mat))+10])
set(gca,'FontSize',8)

% ESPVR 
subplot(3,2,3)
hold on 
h1 = plot(ESV_LV_vec,ESP_LV_vec,'ro','Markersize',6,'MarkerFaceColor','r');
h2 = plot(ESV_RV_vec,ESP_RV_vec,'bo','Markersize',6,'MarkerFaceColor','b');  
l_norm_h = plot(ESV_LV_vec(base_loc),ESP_LV_vec(base_loc),'ko','Markersize',10,'linewidth',2);
r_norm_h = plot(ESV_RV_vec(base_loc),ESP_RV_vec(base_loc),'ko','Markersize',10,'linewidth',2);
legend([h1 h2 l_norm_h],'LV','RV','Baseline','Location','northeast')
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')    
set(gca,'FontSize',8)

% EDPVR 
subplot(3,2,4)
hold on 
plot(V_EDPVR,P_LV_EDPVR_h,'color',gray,'linewidth',2)
plot(V_EDPVR,P_RV_EDPVR_h,'color',gray,'linewidth',2)
h1 = plot(EDV_LV_vec,EDP_LV_vec,'ro','Markersize',6,'MarkerFaceColor','r');
h2 = plot(EDV_RV_vec,EDP_RV_vec,'bo','Markersize',6,'MarkerFaceColor','b');
l_norm_h = plot(EDV_LV_vec(base_loc),EDP_LV_vec(base_loc),'ko','Markersize',10,'linewidth',2);
r_norm_h = plot(EDV_RV_vec(base_loc),EDP_RV_vec(base_loc),'ko','Markersize',10,'linewidth',2);
legend([h1 h2 l_norm_h],'LV','RV','Baseline','Location','northwest')
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
xlim([50 250])
ylim([0 50])
set(gca,'FontSize',8)

% Frank Starling, LV 
subplot(3,2,5)
hold on
plot(EDP_LV_vec, SV_LV_vec,'o','color',green','MarkerSize',6,'MarkerFaceColor',green)
h_norm_h = plot(EDP_LV_vec(base_loc),SV_LV_vec(base_loc),'ko','Markersize',10,'linewidth',2);
legend([h_norm_h],'Baseline','Location','southeast')
xlabel('EDP (mmHg)')
ylabel('SV (mL)')
set(gca,'FontSize',8)

% Frank Starling, RV
subplot(3,2,6)
hold on
plot(EDP_RV_vec, SV_RV_vec,'o','color',green','MarkerSize',6,'MarkerFaceColor',green)
h_norm_h = plot(EDP_RV_vec(base_loc),SV_RV_vec(base_loc),'ko','Markersize',10,'linewidth',2);
legend([h_norm_h],'Baseline','Location','southeast')
xlabel('EDP (mmHg)')
ylabel('SV (mL)')
set(gca,'FontSize',8)

%% Print figures

if printoutfigs_on == 1
    if ~exist('Figures', 'dir')
        mkdir('Figures')
    end
    if ~exist(strcat('Figures/',run_experiment), 'dir')
        mkdir(strcat('Figures/',run_experiment))
    end
    print(hfig1,'-dpng',strcat('Figures/',run_experiment,'/F3a_PVloops.png'))
    print(hfig2,'-dpng',strcat('Figures/',run_experiment,'/F3b_Volumes.png'))
    print(hfig3,'-dpng',strcat('Figures/',run_experiment,'/F3c_HighPressures.png'))
    print(hfig4,'-dpng',strcat('Figures/',run_experiment,'/F3d_LowPressures.png'))
    print(hfig5,'-dpng',strcat('Figures/',run_experiment,'/F_supp_Flows.png'))
    print(hfig6,'-dpng',strcat('Figures/',run_experiment,'/F_supp_Activation.png'))
    print(hfig7,'-dpng',strcat('Figures/',run_experiment,'/F_supp_TBV.png'))
    print(hfig8,'-dpng',strcat('Figures/',run_experiment,'/F4a_LV_progression.png'))
    print(hfig9,'-dpng',strcat('Figures/',run_experiment,'/F4b_RV_progression.png'))
    print(hfig10,'-dpng',strcat('Figures/',run_experiment,'/F4c_EDPVR.png'))
    print(hfig11,'-dpng',strcat('Figures/',run_experiment,'/F4d_ESPVR.png'))
    print(hfig12,'-dpng',strcat('Figures/',run_experiment,'/F4e_LV_progression.png'))
    print(hfig13,'-dpng',strcat('Figures/',run_experiment,'/F4f_RV_progression.png'))
    print(hfig14,'-dpng',strcat('Figures/',run_experiment,'/F4e_SEP_curvature.png'))
    print(hfig15,'-dpng',strcat('Figures/','ExVivo','/F2_g_LV_EDPVR.png'))
    print(hfig16,'-dpng',strcat('Figures/','ExVivo','/F2_h_RV_EDPVR.png'))
    print(hfig17,'-dpng',strcat('Figures/',run_experiment,'/Figure_3.png'))
    print(hfig18,'-dpng',strcat('Figures/',run_experiment,'/Figure_4.png'))


    print(hfig1,'-depsc2',strcat('Figures/',run_experiment,'/F3a_PVloops.eps'))
    print(hfig2,'-depsc2',strcat('Figures/',run_experiment,'/F3b_Volumes.eps'))
    print(hfig3,'-depsc2',strcat('Figures/',run_experiment,'/F3c_HighPressures.eps'))
    print(hfig4,'-depsc2',strcat('Figures/',run_experiment,'/F3d_LowPressures.eps'))
    print(hfig5,'-depsc2',strcat('Figures/',run_experiment,'/F_supp_Flows.eps'))
    print(hfig6,'-depsc2',strcat('Figures/',run_experiment,'/F_supp_Activation.eps'))
    print(hfig7,'-depsc2',strcat('Figures/',run_experiment,'/F_supp_TBV.eps'))
    print(hfig8,'-depsc2',strcat('Figures/',run_experiment,'/F4a_LV_progression.eps'))
    print(hfig9,'-depsc2',strcat('Figures/',run_experiment,'/F4b_RV_progression.eps'))
    print(hfig10,'-depsc2',strcat('Figures/',run_experiment,'/F4c_EDPVR.eps'))
    print(hfig11,'-depsc2',strcat('Figures/',run_experiment,'/F4d_ESPVR.eps'))
    print(hfig12,'-depsc2',strcat('Figures/',run_experiment,'/F4e_LV_progression.eps'))
    print(hfig13,'-depsc2',strcat('Figures/',run_experiment,'/F4f_RV_progression.eps'))
    print(hfig14,'-depsc2',strcat('Figures/',run_experiment,'/F4e_SEP_curvature.eps'))
    print(hfig15,'-depsc2',strcat('Figures/','ExVivo','/F2_g_LV_EDPVR.eps'))
    print(hfig16,'-depsc2',strcat('Figures/','ExVivo','/F2_h_RV_EDPVR.eps'))
    print(hfig17,'-depsc2',strcat('Figures/',run_experiment,'/Figure_3.eps'))
    set(hfig18,'renderer','painters')
    print(hfig18,'-depsc2',strcat('Figures/',run_experiment,'/Figure_4.eps'))

end

function [] = plot_invivo_VL_healthymodsev(loading,data)

%{ 

Plot the volume loading experiment for the in vivo version of the model
for the healthy case and the moderate and severe cases for a particular
dysfunction. 

Inputs: 
loading         - input structure from volumeloading.m containing 
                      vectors for clinical metrics for each case 

data            - input data structure with data and global parameters 
    
Outputs: (figures) 

%} 

%% Flags

printoutfigs_on = data.printoutfigs_on;

%% Unpack loading structure 

run_experiment = loading.run_experiment; 
loading_h = loading.healthy; 
loading_m = loading.moderate; 
loading_s = loading.severe; 

%% Unpack healthy loading structure 

time_h = loading_h.time; 

T = loading_h.T; 

base_loc_h = loading_h.base_loc; 

EDV_LV_vec_h = loading_h.EDV_LV_vec;    EDV_RV_vec_h = loading_h.EDV_RV_vec; 
EDP_LV_vec_h = loading_h.EDP_LV_vec;    EDP_RV_vec_h = loading_h.EDP_RV_vec; 
ESV_LV_vec_h = loading_h.ESV_LV_vec;    ESV_RV_vec_h = loading_h.ESV_RV_vec; 
ESP_LV_vec_h = loading_h.ESP_LV_vec;    ESP_RV_vec_h = loading_h.ESP_RV_vec; 

SV_LV_vec_h = loading_h.SV_LV_vec;
SV_RV_vec_h = loading_h.SV_RV_vec; 

V_LV_mat_h = loading_h.V_LV_mat;    V_RV_mat_h = loading_h.V_RV_mat;
P_LV_mat_h = loading_h.P_LV_mat;    P_RV_mat_h = loading_h.P_RV_mat;

V_LV_h = V_LV_mat_h(:,base_loc_h);  V_RV_h = V_RV_mat_h(:,base_loc_h);
P_LV_h = P_LV_mat_h(:,base_loc_h);  P_RV_h = P_RV_mat_h(:,base_loc_h);

CO_LV_vec_h = loading_h.CO_LV_vec;  CO_RV_vec_h = loading_h.CO_RV_vec; 
CO_LV_h = CO_LV_vec_h(base_loc_h);  CO_RV_h = CO_RV_vec_h(base_loc_h); 

CP_LV_vec_h = loading_h.CP_LV_vec;  CP_RV_vec_h = loading_h.CP_RV_vec; 
CP_LV_h = CP_LV_vec_h(base_loc_h);  CP_RV_h = CP_RV_vec_h(base_loc_h);

EF_LV_vec_h = loading_h.EF_LV_vec;  EF_RV_h_vec = loading_h.EF_RV_vec;
EF_LV_h = EF_LV_vec_h(base_loc_h);  EF_RV_h = EF_RV_h_vec(base_loc_h);

Cm_SEP_mat_h = loading_h.Cm_SEP_mat; 

Cm_SEP_vec_h = Cm_SEP_mat_h(:,base_loc_h);

Cm_LV_mat_h = loading_h.Cm_LV_mat; 
Cm_LV_vec_h = Cm_LV_mat_h(:,base_loc_h);
Cm_RV_mat_h = loading_h.Cm_RV_mat; 
Cm_RV_vec_h = Cm_RV_mat_h(:,base_loc_h);

P_peri_mat_h = loading_h.P_peri_mat; 
P_peri_h = P_peri_mat_h(:,base_loc_h);

Q_m_valve_mat_h = loading_h.Q_m_valve_mat; 
Q_m_valve = Q_m_valve_mat_h(:,base_loc_h);

Q_a_valve_mat_h = loading_h.Q_a_valve_mat; 
Q_a_valve = Q_a_valve_mat_h(:,base_loc_h);

% Find ejection
i_AVO = find(diff(Q_a_valve) > 0,1,'first'); 
i_AVC = find(diff(Q_a_valve) < 0,1,'last'); 

T_AVO = time_h(i_AVO); 
T_AVC = time_h(i_AVC) - T; 

% Find filling
i_MVO = find(diff(Q_m_valve) > 0,1,'first'); 
i_MVC = find(diff(Q_m_valve) < 0,1,'last'); 

T_MVO = time_h(i_MVO); 
T_MVC = time_h(i_MVC) - T; 


%% Unpack moderate dysfunction loading structure 

time_m = loading_m.time; 

base_loc_m = loading_m.base_loc; 

EDV_LV_vec_m = loading_m.EDV_LV_vec;    EDV_RV_vec_m = loading_m.EDV_RV_vec; 
EDP_LV_vec_m = loading_m.EDP_LV_vec;    EDP_RV_vec_m = loading_m.EDP_RV_vec; 
ESV_LV_vec_m = loading_m.ESV_LV_vec;    ESV_RV_vec_m = loading_m.ESV_RV_vec; 
ESP_LV_vec_m = loading_m.ESP_LV_vec;    ESP_RV_vec_m = loading_m.ESP_RV_vec; 

SV_LV_vec_m = loading_m.SV_LV_vec;      SV_RV_vec_m = loading_m.SV_RV_vec; 

V_LV_mat_m = loading_m.V_LV_mat;    V_RV_mat_m = loading_m.V_RV_mat;
P_LV_mat_m = loading_m.P_LV_mat;    P_RV_mat_m = loading_m.P_RV_mat;

V_LV_m = V_LV_mat_m(:,base_loc_m);  V_RV_m = V_RV_mat_m(:,base_loc_m);
P_LV_m = P_LV_mat_m(:,base_loc_m);  P_RV_m = P_RV_mat_m(:,base_loc_m);

CO_LV_vec_m = loading_m.CO_LV_vec;  CO_RV_vec_m = loading_m.CO_RV_vec; 
CO_LV_m = CO_LV_vec_m(base_loc_m);  CO_RV_m = CO_RV_vec_m(base_loc_m); 

CP_LV_vec_m = loading_m.CP_LV_vec;  CP_RV_vec_m = loading_m.CP_RV_vec;
CP_LV_m = CP_LV_vec_m(base_loc_m);  CP_RV_m = CP_RV_vec_m(base_loc_m);

EF_LV_vec_m = loading_m.EF_LV_vec;  EF_RV_m_vec = loading_m.EF_RV_vec;
EF_LV_m = EF_LV_vec_m(base_loc_m);  EF_RV_m = EF_RV_m_vec(base_loc_m);

% Cm_SEP_mat_m = loading_m.Cm_SEP_mat; 
% Cm_SEP_vec_m = Cm_SEP_mat_m(:,base_loc_m);

P_peri_mat_m = loading_m.P_peri_mat; 
P_peri_m = P_peri_mat_m(:,base_loc_m);

%% Unpack severe dysfunction loading structure 

time_s = loading_s.time; 

base_loc_s = loading_s.base_loc; 

EDV_LV_vec_s = loading_s.EDV_LV_vec;    EDV_RV_vec_s = loading_s.EDV_RV_vec; 
EDP_LV_vec_s = loading_s.EDP_LV_vec;    EDP_RV_vec_s = loading_s.EDP_RV_vec; 
ESV_LV_vec_s = loading_s.ESV_LV_vec;    ESV_RV_vec_s = loading_s.ESV_RV_vec; 
ESP_LV_vec_s = loading_s.ESP_LV_vec;    ESP_RV_vec_s = loading_s.ESP_RV_vec; 

SV_LV_vec_s = loading_s.SV_LV_vec;      SV_RV_vec_s = loading_s.SV_RV_vec; 

V_LV_mat_s = loading_s.V_LV_mat;        V_RV_mat_s = loading_s.V_RV_mat;
P_LV_mat_s = loading_s.P_LV_mat;        P_RV_mat_s = loading_s.P_RV_mat;

V_LV_s = V_LV_mat_s(:,base_loc_s);      V_RV_s = V_RV_mat_s(:,base_loc_s);
P_LV_s = P_LV_mat_s(:,base_loc_s);      P_RV_s = P_RV_mat_s(:,base_loc_s);

CO_LV_vec_s = loading_s.CO_LV_vec;      CO_RV_vec_s = loading_s.CO_RV_vec; 
CO_LV_s = CO_LV_vec_s(base_loc_s);      CO_RV_s = CO_RV_vec_s(base_loc_s);

CP_LV_vec_s = loading_s.CP_LV_vec;      CP_RV_vec_s = loading_s.CP_RV_vec;
CP_LV_s = CP_LV_vec_s(base_loc_s);      CP_RV_s = CP_RV_vec_s(base_loc_s);

EF_LV_vec_s = loading_s.EF_LV_vec;      EF_RV_s_vec = loading_s.EF_RV_vec;
EF_LV_s = EF_LV_vec_s(base_loc_s);      EF_RV_s = EF_RV_s_vec(base_loc_s);

Cm_SEP_mat_s = loading_s.Cm_SEP_mat; 
Cm_SEP_vec_s = Cm_SEP_mat_s(:,base_loc_s);
Cm_LV_mat_s = loading_s.Cm_LV_mat; 
Cm_LV_vec_s = Cm_LV_mat_s(:,base_loc_s);
Cm_RV_mat_s = loading_s.Cm_RV_mat; 
Cm_RV_vec_s = Cm_RV_mat_s(:,base_loc_s);

P_peri_mat_s = loading_s.P_peri_mat; 
P_peri_s = P_peri_mat_s(:,base_loc_s);

CP_LV = [CP_LV_h; CP_LV_m; CP_LV_s];
CP_RV = [CP_RV_h; CP_RV_m; CP_RV_s];

%% Make EDPVR Klotz curves

V_EDPVR = [10:200];

% [EDV_h,EDP_h] = getEDESvals(V_LV_h,V_RV_h,P_LV_h,P_RV_h); 

EDV_LV_h = EDV_LV_vec_h(base_loc_h);    EDV_RV_h = EDV_RV_vec_h(base_loc_h);
EDP_LV_h = EDP_LV_vec_h(base_loc_h);    EDP_RV_h = EDP_RV_vec_h(base_loc_h);
EDV_h = [EDV_LV_h, EDV_RV_h];
EDP_h = [EDP_LV_h, EDP_RV_h];

[P_LV_EDPVR_h,P_RV_EDPVR_h] = makeKlotzcurve(EDV_h,EDP_h,V_EDPVR);

%% ESPVR curves

xfit = linspace(30,100); 
    
p_LV_h = polyfit(ESV_LV_vec_h(2:end),ESP_LV_vec_h(2:end),1);
p_RV_h = polyfit(ESV_RV_vec_h(2:end),ESP_RV_vec_h(2:end),1);
yfit_LV_h = p_LV_h(1)*xfit + p_LV_h(2); 
yfit_RV_h = p_RV_h(1)*xfit + p_RV_h(2); 

p_LV_m = polyfit(ESV_LV_vec_m(2:end),ESP_LV_vec_m(2:end),1);
p_RV_m = polyfit(ESV_RV_vec_m(2:end),ESP_RV_vec_m(2:end),1);
yfit_LV_m = p_LV_m(1)*xfit + p_LV_m(2); 
yfit_RV_m = p_RV_m(1)*xfit + p_RV_m(2); 

p_LV_s = polyfit(ESV_LV_vec_s(2:end),ESP_LV_vec_s(2:end),1);
p_RV_s = polyfit(ESV_RV_vec_s(2:end),ESP_RV_vec_s(2:end),1);
yfit_LV_s = p_LV_s(1)*xfit + p_LV_s(2);
yfit_RV_s = p_RV_s(1)*xfit + p_RV_s(2); 
    

%% Plot figures

gray = [.5 .5 .5]; 
green = [0 .5 0]; 

%% Healthy moderate severe PV loops at same Vtot 
% LV 
hfig1 = figure(1); 
clf
hold on 
plot(V_EDPVR,P_LV_EDPVR_h,'color',gray,'linewidth',1)
h1 = plot(V_LV_mat_h(:,base_loc_h), P_LV_mat_h(:,base_loc_h), 'r',   'linewidth',2);
h2 = plot(V_LV_mat_m(:,base_loc_h), P_LV_mat_m(:,base_loc_h), 'r:',  'linewidth',2);
h3 = plot(V_LV_mat_s(:,base_loc_h), P_LV_mat_s(:,base_loc_h), 'r--', 'linewidth',2);
% title(strcat(run_experiment,'-LV-healthy volume'))
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
legend([h1 h2 h3],'H','M','S','location','Northwest')

xlim([0 200])
ylim([0 150])
set(gca,'FontSize',20)

% RV  
hfig2 = figure(2); 
clf
hold on 
plot(V_EDPVR,P_RV_EDPVR_h,'color',gray,'linewidth',1)
h1 = plot(V_RV_mat_h(:,base_loc_h), P_RV_mat_h(:,base_loc_h), 'b',   'linewidth',2);
h2 = plot(V_RV_mat_m(:,base_loc_h), P_RV_mat_m(:,base_loc_h), 'b:',  'linewidth',2);
h3 = plot(V_RV_mat_s(:,base_loc_h), P_RV_mat_s(:,base_loc_h), 'b--', 'linewidth',2);

% title(strcat('RV-',run_experiment))

% title(strcat(run_experiment,'-RV-healthy volume'))
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
legend([h1 h2 h3],'H','M','S','location','Northwest')

xlim([0 200])
ylim([0 60])
set(gca,'FontSize',20)

%% LV PV loops 
hfig3 = figure(3); 
clf
hold on 
plot(V_EDPVR,P_LV_EDPVR_h,'color',gray,'linewidth',1);
h1 = plot(V_LV_h, P_LV_h, 'r',   'linewidth',2);
h2 = plot(V_LV_m, P_LV_m, 'r:',  'linewidth',2);
h3 = plot(V_LV_s, P_LV_s, 'r--', 'linewidth',2);

% title(strcat(run_experiment,'-LV'))
% title('LV')
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
legend([h1 h2 h3],'H','M','S','location','Northwest')

xlim([0 200])
ylim([0 150])
set(gca,'FontSize',20)

%% RV PV loops 
hfig4 = figure(4); 
clf
hold on 
plot(V_EDPVR,P_RV_EDPVR_h,'color',gray,'linewidth',1)
h1 = plot(V_RV_h, P_RV_h, 'b',   'linewidth',2);
h2 = plot(V_RV_m, P_RV_m, 'b:',  'linewidth',2);
h3 = plot(V_RV_s, P_RV_s, 'b--', 'linewidth',2);

% title(strcat(run_experiment,'-RV'))

xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
legend([h1 h2 h3],'H','M','S','location','Northwest')

xlim([0 200])
ylim([0 65])
set(gca,'FontSize',20)

%% ESPVR Plots
hfig5 = figure(5);
clf
hold on 
h1 = plot(ESV_LV_vec_h,ESP_LV_vec_h,'ro','Markersize',6,'MarkerFaceColor','r');
h2 = plot(ESV_RV_vec_h,ESP_RV_vec_h,'bo','Markersize',6,'MarkerFaceColor','b');
h3 = plot(ESV_LV_vec_m,ESP_LV_vec_m,'r^','Markersize',6,'MarkerFaceColor','r');
h4 = plot(ESV_RV_vec_m,ESP_RV_vec_m,'b^','Markersize',6,'MarkerFaceColor','b');
h5 = plot(ESV_LV_vec_s,ESP_LV_vec_s,'rs','Markersize',6,'MarkerFaceColor','r');
h6 = plot(ESV_RV_vec_s,ESP_RV_vec_s,'bs','Markersize',6,'MarkerFaceColor','b');
    
l_norm_h = plot(ESV_LV_vec_h(base_loc_h),ESP_LV_vec_h(base_loc_h),'ko','Markersize',10,'linewidth',2);
r_norm_h = plot(ESV_RV_vec_h(base_loc_h),ESP_RV_vec_h(base_loc_h),'ko','Markersize',10,'linewidth',2);
l_norm_m = plot(ESV_LV_vec_m(base_loc_m),ESP_LV_vec_m(base_loc_m),'k^','Markersize',10,'linewidth',2);
r_norm_m = plot(ESV_RV_vec_m(base_loc_m),ESP_RV_vec_m(base_loc_m),'k^','Markersize',10,'linewidth',2);
l_norm_s = plot(ESV_LV_vec_s(base_loc_s),ESP_LV_vec_s(base_loc_s),'ks','Markersize',10,'linewidth',2);
r_norm_s = plot(ESV_RV_vec_s(base_loc_s),ESP_RV_vec_s(base_loc_s),'ks','Markersize',10,'linewidth',2);

legend([l_norm_h, l_norm_m, l_norm_s],'H','M','S')

xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')    
% title(strcat(run_experiment,'-ESPVR'))   

set(gca,'FontSize',20)
xlim([30 150])
ylim([0 150])

%% EDPVR Plots
hfig6 = figure(6);
clf
hold on 
h1 = plot(EDV_LV_vec_h,EDP_LV_vec_h,'ro','Markersize',6,'MarkerFaceColor','r');
h2 = plot(EDV_RV_vec_h,EDP_RV_vec_h,'bo','Markersize',6,'MarkerFaceColor','b');
h3 = plot(EDV_LV_vec_m,EDP_LV_vec_m,'r^','Markersize',6,'MarkerFaceColor','r');
h4 = plot(EDV_RV_vec_m,EDP_RV_vec_m,'b^','Markersize',6,'MarkerFaceColor','b');
h5 = plot(EDV_LV_vec_s,EDP_LV_vec_s,'rs','Markersize',6,'MarkerFaceColor','r');
h6 = plot(EDV_RV_vec_s,EDP_RV_vec_s,'bs','Markersize',6,'MarkerFaceColor','b');
    
l_norm_h = plot(EDV_LV_vec_h(base_loc_h),EDP_LV_vec_h(base_loc_h),'ko','Markersize',10,'linewidth',2);
r_norm_h = plot(EDV_RV_vec_h(base_loc_h),EDP_RV_vec_h(base_loc_h),'ko','Markersize',10,'linewidth',2);
l_norm_m = plot(EDV_LV_vec_m(base_loc_m),EDP_LV_vec_m(base_loc_m),'k^','Markersize',10,'linewidth',2);
r_norm_m = plot(EDV_RV_vec_m(base_loc_m),EDP_RV_vec_m(base_loc_m),'k^','Markersize',10,'linewidth',2);
l_norm_s = plot(EDV_LV_vec_s(base_loc_s),EDP_LV_vec_s(base_loc_s),'ks','Markersize',10,'linewidth',2);
r_norm_s = plot(EDV_RV_vec_s(base_loc_s),EDP_RV_vec_s(base_loc_s),'ks','Markersize',10,'linewidth',2);

legend([l_norm_h, l_norm_m, l_norm_s],'H','M','S','Location','northwest')

xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
% title(strcat(run_experiment,'-EDPVR'))
    
xlim([50 250])
ylim([0 50])

set(gca,'FontSize',20)

    
%% Frank Starling
% LV 
hfig7 = figure(7);
clf
hold on

plot(EDP_LV_vec_h, SV_LV_vec_h,'o','color',green','MarkerSize',6,'MarkerFaceColor',green)
plot(EDP_LV_vec_m, SV_LV_vec_m,'^','color',green','MarkerSize',6,'MarkerFaceColor',green)
plot(EDP_LV_vec_s, SV_LV_vec_s,'s','color',green','MarkerSize',6,'MarkerFaceColor',green)
    
h_norm_h = plot(EDP_LV_vec_h(base_loc_h),SV_LV_vec_h(base_loc_h),'ko','Markersize',10,'linewidth',2);
h_norm_m = plot(EDP_LV_vec_m(base_loc_m),SV_LV_vec_m(base_loc_m),'k^','Markersize',10,'linewidth',2);
h_norm_s = plot(EDP_LV_vec_s(base_loc_s),SV_LV_vec_s(base_loc_s),'ks','Markersize',10,'linewidth',2);

legend([h_norm_h, h_norm_m, h_norm_s],'H','M','S','Location','southeast')
 
xlabel('EDP (mmHg)')
ylabel('SV (mL)')
% title(strcat(run_experiment,'- Frank-Starling - LV')) 
set(gca,'FontSize',20)


% RV
hfig8 = figure(8);
clf
hold on

plot(EDP_RV_vec_h, SV_RV_vec_h,'o','color',green','MarkerSize',6,'MarkerFaceColor',green)
plot(EDP_RV_vec_m, SV_RV_vec_m,'^','color',green','MarkerSize',6,'MarkerFaceColor',green)
plot(EDP_RV_vec_s, SV_RV_vec_s,'s','color',green','MarkerSize',6,'MarkerFaceColor',green)
    
h_norm_h = plot(EDP_RV_vec_h(base_loc_h),SV_RV_vec_h(base_loc_h),'ko','Markersize',10,'linewidth',2);
h_norm_m = plot(EDP_RV_vec_m(base_loc_m),SV_RV_vec_m(base_loc_m),'k^','Markersize',10,'linewidth',2);
h_norm_s = plot(EDP_RV_vec_s(base_loc_s),SV_RV_vec_s(base_loc_s),'ks','Markersize',10,'linewidth',2);

legend([h_norm_h, h_norm_m, h_norm_s],'H','M','S','Location','southeast')
    
xlabel('EDP (mmHg)')
ylabel('SV (mL)')
% title(strcat(run_experiment,'- Frank-Starling - RV')) 
set(gca,'FontSize',20)


%% Cardiac power 
hfig10 = figure(10); 
clf
hold on 
h1 = plot([1:length(CP_LV)],CP_LV,'rp-','MarkerSize',10,'MarkerFaceColor','r');
h2 = plot([1:length(CP_RV)],CP_RV,'bp-','MarkerSize',10,'MarkerFaceColor','b');
set(gca,'XTick',[1 2 3])
set(gca,'XTickLabels',{'H','M','S'})
ylabel('Cardiac power output (W)')

x = [0 4 4 0];
y1 = [.85 .85 1.15 1.15];
patch(x,y1,'r','linestyle','none')
alpha(.05)

y2 = [.2 .2 .3 .3];
patch(x,y2,'b','linestyle','none')
alpha(.05)

legend([h1,h2],'LV','RV')

% title(strcat(run_experiment)) 

xlim([0 4])
ylim([0 1.2])
set(gca,'FontSize',20)

%% Septal Curvature 
hfig9 = figure(9);
clf
hold on 
h1 = plot(time_h,Cm_SEP_vec_h,'color', green,'linewidth',2);
% h2 = plot(time_m,Cm_SEP_vec_m, ':k','color', green,'linewidth',2);
h3 = plot(time_s,Cm_SEP_vec_s, '--','color', green,'linewidth',2);
ylim([0 .4])
xlim([0 1])

x = [T_AVO T_AVC T_AVC T_AVO];
y1 = [0 0 .5 .5];
patch(x,y1,gray,'linestyle','none')
 alpha(.1)
text(T_AVO+.02,.04,'Ejection','fontsize',20)

x = [T_MVO T_MVC T_MVC T_MVO];
y1 = [0 0 .5 .5];
patch(x,y1,gray,'linestyle','none')
alpha(.1)
text(T_MVO+.02,.04,'Filling','fontsize',20)


xlabel('Time (s)')
ylabel('Septal Curvature (cm^{-1})')
% title(strcat(run_experiment)) 
legend([h1 h3],'H','S')
    

set(gca,'FontSize',20)


%% Print figures

    if printoutfigs_on == 1
        if ~exist('Figures', 'dir')
            mkdir('Figures')
        end
        if ~exist(strcat('Figures/',run_experiment), 'dir')
            mkdir(strcat('Figures/',run_experiment))
        end

        print(hfig1, '-dpng',strcat('Figures/',run_experiment,'/Fa_PVloops_LV_HTBV.png'))
        print(hfig2, '-dpng',strcat('Figures/',run_experiment,'/Fb_PVloops_RV_HTBV.png'))
        print(hfig3, '-dpng',strcat('Figures/',run_experiment,'/Fc_PVloops_LV.png'))
        print(hfig4, '-dpng',strcat('Figures/',run_experiment,'/Fd_PVloops_RV.png'))
        print(hfig5, '-dpng',strcat('Figures/',run_experiment,'/Fe_ESPVR.png'))
        print(hfig6, '-dpng',strcat('Figures/',run_experiment,'/Ff_EDPVR.png'))
        print(hfig7, '-dpng',strcat('Figures/',run_experiment,'/Fg_FS_LV.png'))
        print(hfig8, '-dpng',strcat('Figures/',run_experiment,'/Fh_FS_RV.png'))
        print(hfig9, '-dpng',strcat('Figures/',run_experiment,'/Fi_Cm_SEP.png'))
        print(hfig10, '-dpng',strcat('Figures/',run_experiment,'/Fj_CP.png'))

        print(hfig1, '-depsc2',strcat('Figures/',run_experiment,'/Fa_PVloops_LV_HTBV.eps'))
        print(hfig2, '-depsc2',strcat('Figures/',run_experiment,'/Fb_PVloops_RV_HTBV.eps'))
        print(hfig3, '-depsc2',strcat('Figures/',run_experiment,'/Fc_PVloops_LV.eps'))
        print(hfig4, '-depsc2',strcat('Figures/',run_experiment,'/Fd_PVloops_RV.eps'))
        print(hfig5, '-depsc2',strcat('Figures/',run_experiment,'/Fe_ESPVR.eps'))
        print(hfig6, '-depsc2',strcat('Figures/',run_experiment,'/Ff_EDPVR.eps'))
        print(hfig7, '-depsc2',strcat('Figures/',run_experiment,'/Fg_FS_LV.eps'))
        print(hfig8, '-depsc2',strcat('Figures/',run_experiment,'/Fh_FS_RV.eps'))
        print(hfig9, '-depsc2',strcat('Figures/',run_experiment,'/Fi_Cm_SEP.eps'))
        print(hfig10, '-depsc2',strcat('Figures/',run_experiment,'/Fj_CP.eps'))        

    end


end
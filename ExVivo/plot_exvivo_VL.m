function [] = plot_exvivo_VL(outputs,data)


%{ 
    This function makes the plots for the ex vivo model. 
    Inputs: 
    outputs     - output structure from model_sol_exvivo.m 
    data        - input data structure with data and global parameters 
%} 

a_EDP_LV = outputs.EDPVRs.a_EDP_LV; 
a_EDV_LV = outputs.EDPVRs.a_EDV_LV; 
a_EDP_RV = outputs.EDPVRs.a_EDP_RV; 
a_EDV_RV = outputs.EDPVRs.a_EDV_RV; 

V_LV_EDPVR  = outputs.EDPVRs.V_LV_EDPVR;
V_RV_EDPVR  = outputs.EDPVRs.V_RV_EDPVR;
P_LV_EDPVR  = outputs.EDPVRs.P_LV_EDPVR;
P_RV_EDPVR  = outputs.EDPVRs.P_RV_EDPVR;

%% Plot
gray = [.5 .5 .5]; 


hfig1 = figure(1); 
clf
hold on 
h1 = plot(V_LV_EDPVR, P_LV_EDPVR,'color',gray,'linewidth',2); 
h2 = plot(a_EDV_LV, a_EDP_LV, 'ro', 'markersize',6,'markerfacecolor','r'); 
legend([h1 h2],'Klotz','Model','location','northwest')
set(gca,'FontSize',20)
% title('ExVivo EDPVR')
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
ylim([0 50])
xlim([100 200])

hfig2 = figure(2); 
clf
hold on 
h1 = plot(V_RV_EDPVR, P_RV_EDPVR,'color',gray,'linewidth',2); 
h2 = plot(a_EDV_RV, a_EDP_RV, 'bo', 'markersize',6,'markerfacecolor','b'); 
legend([h1 h2],'Klotz','Model','location','northwest')
set(gca,'FontSize',20)
% title('ExVivo EDPVR')
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
ylim([0 50])
xlim([100 200])

%% Print Figures

% print(hfig1,'-dpng',strcat('Figures/','/ExVivo','/F2_a.png'))
% print(hfig2,'-dpng',strcat('Figures/','/ExVivo','/F2_b.png'))
% 
% print(hfig1,'-depsc2',strcat('Figures/','/ExVivo','/F2_a.eps'))
% print(hfig2,'-depsc2',strcat('Figures/','/ExVivo','/F2_b.eps'))


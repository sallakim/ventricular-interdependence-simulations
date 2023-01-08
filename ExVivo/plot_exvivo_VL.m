function [] = plot_exvivo_VL(outputs,data)

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


hfig201 = figure(201); 
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

hfig202 = figure(202); 
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

print(hfig201,'-dpng',strcat('Figures/','/ExVivo','/F2_a_EDPVR_LV.png'))
print(hfig202,'-dpng',strcat('Figures/','/ExVivo','/F2_b_EDPVR_RV.png'))

print(hfig201,'-depsc2',strcat('Figures/','/ExVivo','/F2_a_EDPVR_LV.eps'))
print(hfig202,'-depsc2',strcat('Figures/','/ExVivo','/F2_b_EDPVR_RV.eps'))


% figure(202)
% clf
% plot(rout,'bo')
% txt = strcat('J =',num2str(J)); 
% text(.75*max(length(rout)),.75*max(rout),txt,'FontSize',20)
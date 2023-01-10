function [] = plot_exvivo_VL_EDV(outputs,data)

%{ 
    This function makes the plots for the ex vivo model. 
    Inputs: 
    outputs     - output structure from model_sol_exvivo.m 
    data        - input data structure with data and global parameters 
%} 

% Unpack outputs structure 

    eta_EDV_vec  = outputs.eta_EDV_vec;
    eta_Vtot_vec = data.a_eta_Vtot; 

    EDV_LV_mat = outputs.EDV_LV_mat; 
    EDV_RV_mat = outputs.EDV_RV_mat; 
    EDP_LV_mat = outputs.EDP_LV_mat; 
    EDP_RV_mat = outputs.EDP_RV_mat; 

    EDV_LV_1 = EDV_LV_mat(:,1); EDP_LV_1 = EDP_LV_mat(:,1); 
    EDV_LV_2 = EDV_LV_mat(:,2); EDP_LV_2 = EDP_LV_mat(:,2); 
    EDV_LV_3 = EDV_LV_mat(:,3); EDP_LV_3 = EDP_LV_mat(:,3); 
    EDV_LV_4 = EDV_LV_mat(:,4); EDP_LV_4 = EDP_LV_mat(:,4); 

    EDV_RV_1 = EDV_RV_mat(:,1); EDP_RV_1 = EDP_RV_mat(:,1); 
    EDV_RV_2 = EDV_RV_mat(:,2); EDP_RV_2 = EDP_RV_mat(:,2); 
    EDV_RV_3 = EDV_RV_mat(:,3); EDP_RV_3 = EDP_RV_mat(:,3); 
    EDV_RV_4 = EDV_RV_mat(:,4); EDP_RV_4 = EDP_RV_mat(:,4);

    %% Make Klotz curve 

    EDV_data = [data.EDV_LV/1e-6; data.EDV_RV/1e-6]; 
    EDP_data = [data.EDP_LV*7.5; data.EDP_RV*7.5]; 

    i_EDV = find(eta_EDV_vec == 1); 
    i_Vtot = find(eta_Vtot_vec == 1); 

    V_EDPVR = [EDV_LV_mat(1,i_EDV):2:EDV_LV_mat(end,i_EDV)]'; 

    [P_LV_EDPVR,P_RV_EDPVR,V] = makeKlotzcurve(EDV_data,EDP_data,V_EDPVR); 

    %% Plots

    gray = [.5 .5 .5]; 

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Regular volumes  
    hfig1 = figure(1);
    clf
    hold on 
    h0 = plot(V_EDPVR,P_LV_EDPVR,'color',gray,'linewidth',2);
    h1 = plot(EDV_LV_1,EDP_LV_1,'r-x','linewidth',2);
    h2 = plot(EDV_LV_2,EDP_LV_2,'r-*','linewidth',2);
    h3 = plot(EDV_LV_3,EDP_LV_3,'r-^','linewidth',2);
    h4 = plot(EDV_LV_4,EDP_LV_4,'r-s','linewidth',2);
    ylim([0 50])
    xlim([0 400])
    xlabel('Volume (mL)')
    ylabel('Pressure (mmHg)')
    set(gca,'FontSize',20)
%     title('LV')

    legend([h0,h1,h2,h3,h4],...
        'Klotz', ...
        ['EDV = ', num2str(eta_EDV_vec(1)*data.EDV_LV/1e-6)], ... 
        ['EDV = ', num2str(eta_EDV_vec(2)*data.EDV_LV/1e-6)], ...
        ['EDV = ', num2str(eta_EDV_vec(3)*data.EDV_LV/1e-6)], ...
        ['EDV = ', num2str(eta_EDV_vec(4)*data.EDV_LV/1e-6)])
    
    % RV
    hfig2 = figure(2); 
    clf
    hold on 
    h0 = plot(V_EDPVR,P_RV_EDPVR,'color', gray,'linewidth',2);
    h1 = plot(EDV_RV_1,EDP_RV_1,'b-x','linewidth',2);
    h2 = plot(EDV_RV_2,EDP_RV_2,'b-*','linewidth',2);
    h3 = plot(EDV_RV_3,EDP_RV_3,'b-^','linewidth',2);
    h4 = plot(EDV_RV_4,EDP_RV_4,'b-s','linewidth',2);
    ylim([0 50])
    xlim([0 400])
    xlabel('Volume (mL)')
    ylabel('Pressure (mmHg)')
    set(gca,'FontSize',20)
%     title('RV')
    legend([h0,h1,h2,h3,h4],...
        'Klotz', ...
        ['EDV = ', num2str(eta_EDV_vec(1)*data.EDV_RV/1e-6)], ... 
        ['EDV = ', num2str(eta_EDV_vec(2)*data.EDV_RV/1e-6)], ...
        ['EDV = ', num2str(eta_EDV_vec(3)*data.EDV_RV/1e-6)], ...
        ['EDV = ', num2str(eta_EDV_vec(4)*data.EDV_RV/1e-6)])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Normalized volumes 

    EDV_1 = [EDV_LV_1(i_Vtot); EDV_RV_1(i_Vtot)]; 
    EDV_2 = [EDV_LV_2(i_Vtot); EDV_RV_2(i_Vtot)]; 
    EDV_3 = [EDV_LV_3(i_Vtot); EDV_RV_3(i_Vtot)]; 
    EDV_4 = [EDV_LV_4(i_Vtot); EDV_RV_4(i_Vtot)]; 
    
    EDP_1 = [EDP_LV_1(i_Vtot); EDP_RV_1(i_Vtot)]; 
    EDP_2 = [EDP_LV_2(i_Vtot); EDP_RV_2(i_Vtot)];
    EDP_3 = [EDP_LV_3(i_Vtot); EDP_RV_3(i_Vtot)];
    EDP_4 = [EDP_LV_4(i_Vtot); EDP_RV_4(i_Vtot)];

    [P_LV_EDPVR_1,P_RV_EDPVR_1,V_1] = makeKlotzcurve(EDV_1,EDP_1,V_EDPVR); 
    [P_LV_EDPVR_2,P_RV_EDPVR_2,V_2] = makeKlotzcurve(EDV_2,EDP_2,V_EDPVR); 
    [P_LV_EDPVR_3,P_RV_EDPVR_3,V_3] = makeKlotzcurve(EDV_3,EDP_3,V_EDPVR); 
    [P_LV_EDPVR_4,P_RV_EDPVR_4,V_4] = makeKlotzcurve(EDV_4,EDP_4,V_EDPVR); 

    V_0_LV_1  = V_1(1); 
    V_30_LV_1 = V_1(2); 
    V_0_RV_1  = V_1(3); 
    V_30_RV_1 = V_1(4); 

    V_0_LV_2  = V_2(1);
    V_30_LV_2 = V_2(2);
    V_0_RV_2  = V_2(3);
    V_30_RV_2 = V_2(4);

    V_0_LV_3  = V_3(1);
    V_30_LV_3 = V_3(2);
    V_0_RV_3  = V_3(3);
    V_30_RV_3 = V_3(4);

    V_0_LV_4  = V_4(1);
    V_30_LV_4 = V_4(2);
    V_0_RV_4  = V_4(3);
    V_30_RV_4 = V_4(4);
    
    V_EDPVR_LV_1 = (V_EDPVR - V_0_LV_1) / (V_30_LV_1 - V_0_LV_1);
    V_EDPVR_LV_2 = (V_EDPVR - V_0_LV_2) / (V_30_LV_2 - V_0_LV_2);
    V_EDPVR_LV_3 = (V_EDPVR - V_0_LV_3) / (V_30_LV_3 - V_0_LV_3);
    V_EDPVR_LV_4 = (V_EDPVR - V_0_LV_4) / (V_30_LV_4 - V_0_LV_4);

    V_EDPVR_RV_1 = (V_EDPVR - V_0_RV_1) / (V_30_RV_1 - V_0_RV_1);
    V_EDPVR_RV_2 = (V_EDPVR - V_0_RV_2) / (V_30_RV_2 - V_0_RV_2);
    V_EDPVR_RV_3 = (V_EDPVR - V_0_RV_3) / (V_30_RV_3 - V_0_RV_3);
    V_EDPVR_RV_4 = (V_EDPVR - V_0_RV_4) / (V_30_RV_4 - V_0_RV_4);
   
    hfig3 = figure(3);
    clf
    hold on 
    h1 = plot(V_EDPVR_LV_1,P_LV_EDPVR_1,'r-x','linewidth',2);
    h2 = plot(V_EDPVR_LV_2,P_LV_EDPVR_2,'r-*','linewidth',2);
    h3 = plot(V_EDPVR_LV_3,P_LV_EDPVR_3,'r-^','linewidth',2);
    h4 = plot(V_EDPVR_LV_4,P_LV_EDPVR_4,'r-s','linewidth',2);
    ylim([0 50])
    xlim([0 1.2])
    xlabel('Volume (mL)')
    ylabel('Pressure (mmHg)')
    set(gca,'FontSize',20)
%     title('LV')

    legend([h1,h2,h3,h4],...
        ['EDV = ', num2str(eta_EDV_vec(1)*data.EDV_LV/1e-6)], ... 
        ['EDV = ', num2str(eta_EDV_vec(2)*data.EDV_LV/1e-6)], ...
        ['EDV = ', num2str(eta_EDV_vec(3)*data.EDV_LV/1e-6)], ...
        ['EDV = ', num2str(eta_EDV_vec(4)*data.EDV_LV/1e-6)], ...
        'location','northwest')
    
    hfig4 = figure(4); 
    clf
    hold on 
    h1 = plot(V_EDPVR_RV_1,P_RV_EDPVR_1,'b-x','linewidth',2);
    h2 = plot(V_EDPVR_RV_2,P_RV_EDPVR_2,'b-*','linewidth',2);
    h3 = plot(V_EDPVR_RV_3,P_RV_EDPVR_3,'b-^','linewidth',2);
    h4 = plot(V_EDPVR_RV_4,P_RV_EDPVR_4,'b-s','linewidth',2);
    ylim([0 50])
    xlim([0 1.2])
    xlabel('Volume (mL)')
    ylabel('Pressure (mmHg)')
    set(gca,'FontSize',20)
%     title('RV')
    legend([h1,h2,h3,h4],...
        ['EDV = ', num2str(eta_EDV_vec(1)*data.EDV_RV/1e-6)], ... 
        ['EDV = ', num2str(eta_EDV_vec(2)*data.EDV_RV/1e-6)], ...
        ['EDV = ', num2str(eta_EDV_vec(3)*data.EDV_RV/1e-6)], ...
        ['EDV = ', num2str(eta_EDV_vec(4)*data.EDV_RV/1e-6)], ...
        'location','northwest')




    %% Print figures 

    printoutfigs_on = data.printoutfigs_on; 

    if printoutfigs_on == 1
        if ~exist('Figures', 'dir')
            mkdir('Figures')
        end
        if ~exist('Figures/ExVivo', 'dir')
            mkdir('Figures/ExVivo')
        end

        print(hfig1,'-dpng','Figures/ExVivo/F3_EDPVR_LV_diffEDV.png')
        print(hfig2,'-dpng','Figures/ExVivo/F4_EDPVR_RV_diffEDV.png')
        print(hfig3,'-dpng','Figures/ExVivo/F5_EDPVR_LV_diffEDV_norm.png')
        print(hfig4,'-dpng','Figures/ExVivo/F6_EDPVR_RV_diffEDV_norm.png')

        print(hfig1,'-depsc2','Figures/ExVivo/F3_EDPVR_LV_diffEDV.eps')
        print(hfig2,'-depsc2','Figures/ExVivo/F4_EDPVR_RV_diffEDV.eps')
        print(hfig3,'-depsc2','Figures/ExVivo/F5_EDPVR_LV_diffEDV_norm.eps')
        print(hfig4,'-depsc2','Figures/ExVivo/F6_EDPVR_RV_diffEDV_norm.eps')

    end         

end 


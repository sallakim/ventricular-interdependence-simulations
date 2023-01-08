function [EDV, EDP, ESV, ESP] = getEDESvals(V_LV,V_RV,P_LV,P_RV) % only load in structure instead of these specific volumes a
% and unpack structure from model sol 

% Only look for volumes before tricuspid opens / mitral for edv/p
% In model sol - find where the valves open and close and save those times
% to use here. 
%  

% Calculated model predicted EDP, ESP, EDV, ESV 

% Find max and min volumes for LV and RV 
V_LV_min = min(round(V_LV));    V_LV_max = max(round(V_LV)); 
V_RV_min = min(round(V_RV));    V_RV_max = max(round(V_RV));

% EDV: max volume               % ESV: min volume 
EDV_LV = max(V_LV);             ESV_LV = min(V_LV);
EDV_RV = max(V_RV);             ESV_RV = min(V_RV); 

% EDP: min pressure at max volume 
l_EDP_LV = find(round(V_LV) == V_LV_max); 
EDP_LV = min(P_LV(l_EDP_LV));

% dV_RV = diff(V_LV_d(l_EDP_LV)); 
% ll = find(abs(dV_RV) <= 1e-3,1); 
% EDP_LV = P_LV(l_EDP_LV(ll));
% EDV_LV = V_LV_d(l_EDP_LV(ll));


l_EDP_RV = find(round(V_RV) == V_RV_max); 
EDP_RV = min(P_RV(l_EDP_RV));

% dV_RV = diff(V_RV_d(l_EDP_RV)); 
% ll = find(abs(dV_RV) <= 1e-4,1); 
% EDP_RV = P_RV(l_EDP_RV(ll));
% EDV_RV = V_RV_d(l_EDP_RV(ll));

% ESP: max pressure at min volume 
l_ESP_LV = find(round(V_LV) == V_LV_min);   

% dV_LV = diff(V_LV(l_ESP_LV)); 
% ll = find(abs(dV_LV) <= 1e-3,1); 
% ESP_LV = P_LV(l_ESP_LV(ll+1));
% ESV_LV = V_LV(l_ESP_LV(ll+1));


ESP_LV = max(P_LV(l_ESP_LV));
    
l_ESP_RV = find(round(V_RV) == V_RV_min);
ESP_RV = max(P_RV(l_ESP_RV));

% dV_RV = diff(V_RV(l_ESP_RV)); 
% ll = find(abs(dV_RV) <= 1e-3,1); 
% ESP_RV = P_RV(l_ESP_RV(ll+1));
% ESV_RV = V_RV(l_ESP_RV(ll+1));


% 
% t = [0:length(V_LV)-1]*.001; 
% figure(101)
% clf
% plot(t,V_LV,'k',t(l_ESP_LV),V_LV(l_ESP_LV),'ro')
% yyaxis right 
% plot(t,P_LV,'k',t(l_ESP_LV),P_LV(l_ESP_LV),'ro')


i_ED_LV = find_corner(V_LV,P_LV,'lowerright'); 
i_ED_RV = find_corner(V_RV,P_RV,'lowerright'); 
i_ES_LV = find_corner(V_LV,P_LV,'upperleft'); 
i_ES_RV = find_corner(V_RV,P_RV,'upperleft'); 

EDV_LV = V_LV(i_ED_LV); 
EDV_RV = V_RV(i_ED_RV); 
EDP_LV = P_LV(i_ED_LV); 
EDP_RV = P_RV(i_ED_RV); 

ESV_LV = V_LV(i_ES_LV); 
ESV_RV = V_RV(i_ES_RV); 
ESP_LV = P_LV(i_ES_LV); 
ESP_RV = P_RV(i_ES_RV); 



%% Ouputs 

EDV = [EDV_LV; EDV_RV];
EDP = [EDP_LV; EDP_RV];
ESV = [ESV_LV; ESV_RV];
ESP = [ESP_LV; ESP_RV];

end
    
    

     
    
    
     


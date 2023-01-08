function [outputs,rout,J] = model_sol(adjpars,data)

%% Initialization  

plot_switch = 0; % 0 = off, 1 = on 

% "undo" log from parameters.m
adjpars = exp(adjpars); 

tspan = data.tspan;  
dt    = data.dt; 

HR = data.HR; 
T = 60/HR; 

ODE_TOL = data.gpars.ODE_TOL;

fixpars = data.fixpars;

%% Get initial conditions

if isempty(data.eta_Vtot)
    eta_Vtot = 1; 
else 
    eta_Vtot = data.eta_Vtot; 
end 

init = initialconditions(adjpars,data,eta_Vtot);

%% Set mass matrix M for DAE 
M = speye(length(init));
M(1,1) = 0;
M(2,2) = 0;
M(3,3) = 0;
M(4,4) = 0; 
% m = 15; % Number of states 
% dd = ones(m,1); 
% dd(1:4) = 0; 
% M = spdiags(dd,0,m,m); 
% opts = odeset('Mass',M,'RelTol',ODE_TOL,'AbsTol',ODE_TOL);

% try

% %% Solve model 
% 
% % Use a while loop to allow model to converge to steady-state 
% ndone = 0; 
% while ndone == 0 
% 
%     opts = odeset('Mass',M,'RelTol',ODE_TOL,'AbsTol',ODE_TOL,'MaxStep',1e-1);
%     sol  = ode15s(@model_exvivo,[tspan(1) tspan(end)],init,opts,adjpars,data);
%     
%     % Plot intermediate states if plot_switch is "on"
%     if plot_switch == 1 
%         % Displacements states 1-4
%         figure(111)
%         clf
%         plot(sol.x,sol.y(1:4,:))
%         legend('1','2','3','4')
%         
%         % Sarcomere lengths states 5-7
%         figure(112)
%         clf
%         plot(sol.x,sol.y(5:7,:))
%         legend('5','6','7')
%         
%         % Volumes states 8-11
%         figure(113)
%         clf
%         plot(sol.x,sol.y(8:9,:))
%         legend('8','9')
%         
% 
%     end 
% 
%     if sol.x(end) ~= tspan(end) 
%         % Check to see if the model solved to then of tspan. If not, set the
%         % initial conditions for the next loop at the start of the previous
%         % beat and solve for 10 beats 
%         t = sol.x(1):dt:sol.x(end); 
%         beats = mod(t,T); 
%         x = find(round(beats,3) == 0);
%         y = find(t(x) <= t(end)); 
%         tspan = tspan(1):dt:tspan(x(y(end)));
%     
%         sols = deval(sol,tspan);
%         init  = sols(:,x(y(end-1))); 
%     
%         tspan = tspan(x(y(end-1))):dt:tspan(x(y(end-1))) + 10*T; 
%          
%     else 
%         % If the model has successfully solved at least 10 beats, then we can
%         % assess whether the model has reached steady state 
%         
%         % Extract sarcomere lengths and systemic arterial pressure (Psa)     
%         sols = deval(sol,tspan);
%             
%         Lsc_lv  = sols(5,:) * 1e6;
%         Lsc_sep = sols(6,:) * 1e6; 
%         Lsc_rv  = sols(7,:) * 1e6; 
%         
%         % Plot extracted quantities if plot_switch is "on" 
%         if plot_switch == 1
%             figure(101)
%             clf
%             hold on 
%             plot(tspan,Lsc_lv) 
%             
%             figure(102)
%             clf
%             hold on 
%             plot(tspan,Lsc_sep) 
%             
%             figure(103)
%             clf
%             hold on 
%             plot(tspan,Lsc_rv) 
%             
%         end
%         
%         % Find the last 5 beats of the simulation 
%         xx = find(tspan >= tspan(end) - 5*T); 
%         
%         % Set a peak threshold as half of the amplitude of the last 5 beats 
%         threshold_lv  = (max(Lsc_lv(xx))  - min(Lsc_lv(xx)))/2;
%         threshold_sep = (max(Lsc_sep(xx)) - min(Lsc_sep(xx)))/2;
%         threshold_rv  = (max(Lsc_rv(xx))  - min(Lsc_rv(xx)))/2;
%         
%         % Determine the length of half of the cardiac cycle 
%         half_per = round((T/2)/dt); 
%         
%         % Find peaks for the last 5 beats 
%         [pks_Lsc_lv, loc_pks_Lsc_lv]  = findpeaks(...
%             Lsc_lv,'MinPeakDistance',half_per,'MinPeakProminence',threshold_lv); 
%         [pks_Lsc_sep,loc_pks_Lsc_sep] = findpeaks(...
%             Lsc_sep,'MinPeakDistance',half_per,'MinPeakProminence',threshold_sep); 
%         [pks_Lsc_rv, loc_pks_Lsc_rv]  = findpeaks(...
%             Lsc_rv,'MinPeakDistance',half_per,'MinPeakProminence',threshold_rv); 
%         
%         % Exclude the last peak (so there are 4 peaks)
%         pks_Lsc_lv  = pks_Lsc_lv(end-5:end-1); 
%         pks_Lsc_sep = pks_Lsc_sep(end-5:end-1); 
%         pks_Lsc_rv  = pks_Lsc_rv(end-5:end-1);
%         
%         % Find the locations of the peaks 
%         loc_pks_Lsc_lv  = loc_pks_Lsc_lv(end-5:end-1); 
%         loc_pks_Lsc_sep = loc_pks_Lsc_sep(end-5:end-1); 
%         loc_pks_Lsc_rv  = loc_pks_Lsc_rv(end-5:end-1); 
%         
%         % Find the times where the peaks occur 
%         t_pks_Lsc_lv  = tspan(loc_pks_Lsc_lv);
%         t_pks_Lsc_sep = tspan(loc_pks_Lsc_sep);
%         t_pks_Lsc_rv  = tspan(loc_pks_Lsc_rv);
%         
%         % Create a linear regression through the peaks 
%         pf_Lsc_lv  = polyfit(t_pks_Lsc_lv,pks_Lsc_lv,1); 
%         pf_Lsc_sep = polyfit(t_pks_Lsc_sep,pks_Lsc_sep,1); 
%         pf_Lsc_rv  = polyfit(t_pks_Lsc_rv,pks_Lsc_rv,1); 
%         
%         % Extract the slope of the regression line 
%         slope_Lsc_lv  = pf_Lsc_lv(1);
%         slope_Lsc_sep = pf_Lsc_sep(1);
%         slope_Lsc_rv  = pf_Lsc_rv(1);
%         
%         % Plot regression line through peaks if plot_switch is "on" 
%         if plot_switch == 1
%             % Draw the regression line through the peaks
%             y_Lsc_lv  = polyval(pf_Lsc_lv,tspan); 
%             y_Lsc_sep = polyval(pf_Lsc_sep,tspan); 
%             y_Lsc_rv  = polyval(pf_Lsc_rv,tspan); 
%             
%             % LV sarcomere length 
%             figure(101)
%             hold on 
%             plot(t_pks_Lsc_lv,pks_Lsc_lv,'r*')
%             plot(tspan,y_Lsc_lv,'k')
%             
%             % SEP sarcomere length 
%             figure(102)
%             hold on 
%             plot(t_pks_Lsc_sep,pks_Lsc_sep,'r*')
%             plot(tspan,y_Lsc_sep,'k')
%             
%             % RV sarcomere length 
%             figure(103)
%             hold on 
%             plot(t_pks_Lsc_rv,pks_Lsc_rv,'r*')
%             plot(tspan,y_Lsc_rv,'k')
%             
% 
%         end
%         
%         % If the slope is sufficiently small (i.e. flat), we have reached 
%         % steady state 
%         slope_lim = 1e-3;
%             % Stopping criteria 
%             if  abs(slope_Lsc_lv) < slope_lim && ...
%                     abs(slope_Lsc_sep) < slope_lim && abs(slope_Lsc_rv) < slope_lim
%                 ndone = 1; 
%             end 
%             
%         % If we have not reached steady-state, solve the model for 10 more
%         % beats and reassess convergence 
%         beats = mod(tspan,T); 
%         x = find(round(beats,3) == 0);
%         y = find(tspan(x) <= tspan(end));
%         tspan = tspan(x(y(end-1))):dt:tspan(x(y(end-1))) + 10*T;
%         init  = sols(:,x(y(end-1))); 
%     end 
% end
% 
% % After determining that the model is in steady-state, solve 2 more beats 
% time = [0:dt:2*T]; 
% sol  = ode15s(@model_exvivo,[time(1) time(end)],init,opts,adjpars,data);
% sols = deval(sol,time);
% sols = sols'; 

%% Calculate other time-varying model quantities (pressures, flows, etc.) 

%o = zeros(26,length(time));  
%for i = 1:length(time) 

    [~,o] = model_exvivo([],init,adjpars,data);
    sols = init'; 
%end 

%% Create output structure  

%outputs.time = time; 

% Convert m to cm 
displacements.xm_LV  = sols(:,1) * 1e2; 
displacements.xm_SEP = sols(:,2) * 1e2; 
displacements.xm_RV  = sols(:,3) * 1e2; 
displacements.ym     = sols(:,4) * 1e2; 

% Convert m to um
lengths.Lsc_LV  = sols(:,5) * 1e6;
lengths.Lsc_SEP = sols(:,6) * 1e6; 
lengths.Lsc_RV  = sols(:,7) * 1e6; 

% Convert to m^3 to mL
volumes.V_LV = sols(:,8) * 1e6; 
volumes.V_RV = sols(:,9) * 1e6; 

% Convert kPa to mmHg
pressures.P_LV = o(1,:) * 7.5; 
pressures.P_RV = o(2,:) * 7.5; 

% Convert m^3 to cm^3
wallvolumes.Vm_LV  = o(3,:)  * 1e6; 
wallvolumes.Vm_SEP = o(4,:) * 1e6; 
wallvolumes.Vm_RV  = o(5,:) * 1e6; 

% Convert m^2 to cm^2
areas.Am_LV  = o(6,:) * 1e4; 
areas.Am_SEP = o(7,:) * 1e4; 
areas.Am_RV  = o(8,:) * 1e4; 

% Convert m^(-1) to cm^(-1)
curvatures.Cm_LV  = o(9,:) * 1e-2;
curvatures.Cm_SEP = o(10,:) * 1e-2;
curvatures.Cm_RV  = o(11,:) * 1e-2; 

strains.eps_LV  = o(12,:); 
strains.eps_SEP = o(13,:); 
strains.eps_RV  = o(14,:); 

stresses.passive.sigma_pas_LV  = o(15,:);
stresses.passive.sigma_pas_SEP = o(16,:);
stresses.passive.sigma_pas_RV  = o(17,:);

stresses.active.sigma_act_LV  = o(18,:);
stresses.active.sigma_act_SEP = o(19,:);
stresses.active.sigma_act_RV  = o(20,:);

stresses.total.sigma_LV  = o(21,:);
stresses.total.sigma_SEP = o(22,:);
stresses.total.sigma_RV  = o(23,:);

tensions.Tm_LV = o(24,:);
tensions.Tm_SEP = o(25,:); 
tensions.Tm_RV = o(26,:); 

outputs.volumes       = volumes; 
outputs.pressures     = pressures; 
outputs.displacements = displacements; 
outputs.areas         = areas;
outputs.wallvolumes   = wallvolumes; 
outputs.curvatures    = curvatures; 
outputs.strains       = strains; 
outputs.stresses      = stresses;
outputs.lengths       = lengths; 
outputs.tensions      = tensions; 

%% Sensitivity Analysis 

% % Cardiac output (L min^(-1))
% % range: 5.5 +/- 1
% Q_a_valve = outputs.flows.Q_a_valve;
% CO_model = trapz(time/60,Q_a_valve)/(time(end)/60 - time(1)/60); %SV * HR_end * 1e-3 % L min^(-1)
% CO_data = data.CO*1e6/1000*60; % L / min (input data) 
% if CO_model < 4.5 && CO_model > 6.5
%     rout_a1 = (CO_model - CO_data)/CO_data; 
% else 
%     rout_a1 = .25 * (CO_model - CO_data)/CO_data; 
% end 
% 
% % LV end-diastolic and end-systolic volumes (mL) 
% V_lv = outputs.volumes.V_lv; 
% V_lv_d = max(V_lv); 
% V_lv_s = min(V_lv); 
% V_lv_d_data = data.EDV_lv*1e6; 
% V_lv_s_data = data.ESV_lv*1e6; 
% rout_a2 = [(V_lv_d - V_lv_d_data)/V_lv_d_data; 
%     (V_lv_s - V_lv_s_data)/V_lv_s_data]; 
% 
% % RV end-diastolic and end-systolic volumes (mL) 
% V_rv = outputs.volumes.V_rv; 
% V_rv_d = max(V_rv); 
% V_rv_s = min(V_rv); 
% V_rv_d_data = data.EDV_rv*1e6; 
% V_rv_s_data = data.ESV_rv*1e6; 
% rout_a3 = [(V_rv_d - V_rv_d_data)/V_rv_d_data; 
%     (V_rv_s - V_rv_s_data)/V_rv_s_data];
% 
% % Cardiac power = cardiac power output index / stroke work: 0.5-0.7W/m^2, avg surface area of a human 1.7
% % m^2 --> approx CP = 1 W
% 
% % LV cardiac power (W) 
% % range: 1 +/- 0.2
% P_lv = outputs.pressures.P_lv; 
% CP_model = trapz(P_lv,V_lv) / 7.5 * 1e-3 * HR/60; %mean(P_sa(beat)) / 7.5 * 1e3 * SV * 1e-6 * HR_end/60 % W 
% CP_model = CP_model/2; 
% CP_data = 1; % W
% if CP_model < .8 && CP_model > 1.2
%     rout_a4 = (CP_model - CP_data)/CP_data; 
% else 
%     rout_a4 = .25 * (CP_model - CP_data)/CP_data; 
% end 
% 
% % Diastolic systemic arterial pressure (mmHg)
% % range: 80 +/- 5
% P_sa = outputs.pressures.P_sa; 
% P_sa_d_model = min(P_sa); 
% P_sa_d_data = 80; 
% if P_sa_d_model < 70 && P_sa_d_model > 80
%     rout_p1 = (P_sa_d_model - P_sa_d_data)/P_sa_d_data; 
% else 
%     rout_p1 = .25 * (P_sa_d_model - P_sa_d_data)/P_sa_d_data; 
% end
% 
% % Systolic systemic arterial pressure (mmHg)
% % range: 120 +/- 5
% P_sa_s_model = max(P_sa); 
% P_sa_s_data = 120; 
% if P_sa_s_model < 115 && P_sa_s_model > 125
%     rout_p2 = (P_sa_s_model - P_sa_s_data)/P_sa_s_data; 
% else 
%     rout_p2 = .25 * (P_sa_s_model - P_sa_s_data)/P_sa_s_data;
% end
% 
% % Mean pulmonary arterial pressure (mmHg)
% % range: 14 +/- 3 
% P_pa = outputs.pressures.P_pa; 
% mP_pa_model = trapz(time,P_pa)/(time(end) - time(1));
% mP_pa_data = 14; % mmHg 
% if mP_pa_model < 11 && mP_pa_model > 17
%     rout_p3 = (mP_pa_model - mP_pa_data)/mP_pa_data; 
% else 
%     rout_p3 = .25 * (mP_pa_model - mP_pa_data)/mP_pa_data; 
% end 

% % Right atrial pressure (mmHg)
% % range: 4 +/- 2 
% P_ra = outputs.pressures.P_ra; 
% mP_ra_model = trapz(time,P_ra)/(time(end) - time(1));
% mP_ra_data = 4; % mmHg 
% if mP_ra_model < 2 && mP_ra_model > 6
%     rout_p4 = (mP_ra_model - mP_ra_data)/mP_ra_data; 
% else 
%     rout_p4 = .25 * (mP_ra_model - mP_ra_data)/mP_ra_data; 
% end 

% % Left atrial pressure (mmHg) taken as an estimate of pulmonary capillary
% % wedge pressure (PCWP) 
% % range: 6-12 mmHg 
% P_la = outputs.pressures.P_la; 
% mP_la_model = trapz(time,P_la)/(time(end) - time(1)); 
% mP_la_data = 9; % mmHg 
% if mP_la_model < 6 && mP_la_model > 12
%     rout_p5 = (mP_la_model - mP_la_data)/mP_la_data; 
% else 
%     rout_p5 = .25 * (mP_la_model - mP_la_data)/mP_la_data; 
% % end 
% 
% % Right ventricle systolic pressure (mmHg)
% % range: 25 +/- 5 
% P_rv = outputs.pressures.P_rv; 
% V_rv_min = min(round(V_rv));
% loc = find(round(V_rv) == V_rv_min); 
% P_rv_s_model = max(P_rv(loc)); 
% P_rv_s_data = data.ESP_rv*7.5; % mmHg 
% if P_rv_s_model < 20 && P_rv_s_model > 30 
%     rout_p6 = (P_rv_s_model - P_rv_s_data)/P_rv_s_data; 
% else 
%     rout_p6 = .25 * (P_rv_s_model - P_rv_s_data)/P_rv_s_data; 
% end
% 
% % Right ventricle diastolic pressure
% % range: 5 +/- 2, 0 - 5
% V_rv_max = max(round(V_rv));
% loc = find(round(V_rv) == V_rv_max); 
% P_rv_d_model = min(P_rv(loc)); 
% P_rv_d_data = data.EDP_rv*7.5;% mmHg 
% if P_rv_d_model < 0 && P_rv_d_model > 5
%     rout_p7 = (P_rv_d_model - P_rv_d_data)/P_rv_d_data; 
% else
%     rout_p7 = .25 * (P_rv_d_model - P_rv_d_data)/P_rv_d_data;
% end
% % if P_rv_d_model < 3 && P_rv_d_model > 7 
% %     rout_p7 = (P_rv_d_model - P_rv_d_data)/P_rv_d_data; 
% % else 
% %     rout_p7 = .25 * (P_rv_d_model - P_rv_d_data)/P_rv_d_data; 
% % end
% 
% % LV end-systolic pressure (mmHg)
% % range: 90-150
% V_lv_min = min(round(V_lv));
% loc = find(round(V_lv) == V_lv_min); 
% P_lv_s_model = max(P_lv(loc)); 
% P_lv_s_data = data.ESP_lv*7.5;% mmHg 
% if P_lv_s_model < 110 && P_lv_s_model > 150
%     rout_p8 = (P_lv_s_model - P_lv_s_data)/P_lv_s_data; 
% else 
%     rout_p8 = .25 * (P_lv_s_model - P_lv_s_data)/P_lv_s_data; 
% end
% 
% % LV end-diastolic pressure (mmHg)
% % range: 5-12
% V_lv_max = max(round(V_lv));
% loc = find(round(V_lv) == V_lv_max); 
% P_lv_d_model = min(P_lv(loc)); 
% P_lv_d_data = data.EDP_rv*7.5; % mmHg 
% if P_lv_d_model < 5 && P_lv_d_model > 12 
%     rout_p9 = (P_lv_d_model - P_lv_d_data)/P_lv_d_data; 
% else 
%     rout_p9 = .25 * (P_lv_d_model - P_lv_d_data)/P_lv_d_data; 
% end

% % E/A ratio 
% % range: 1-2
% Q_m_valve = outputs.flows.Q_m_valve;
% [mitral_peaks, loc_mitral_peaks] = findpeaks(Q_m_valve);
% %[mitral_peaks, ~ ] = findpeaks(Q_m_valve(beat));
% E = mitral_peaks(end-1);
% A = mitral_peaks(end);
% EA_model = E/A; 
% EA_data = 1.7; 
% if EA_model < 1 && EA_model > 2
%     rout_f1 = (EA_model - EA_data)/EA_data; 
% else
%     rout_f1 = .25 * (EA_model - EA_data)/EA_data; 
% end
% 
% rout = [rout_a1; rout_a2; rout_a3; rout_a4; 
%         rout_p1; rout_p2; rout_p3; %rout_p4; rout_p5; 
%         rout_p6; rout_p7; rout_p8; rout_p9; 
%         rout_f1]; 
% catch 
%        outputs = [];
%        rout = 100;
% end

rout = []; 

J = rout'*rout;

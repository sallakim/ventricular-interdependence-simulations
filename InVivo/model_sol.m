function [outputs,rout,J] = model_sol(adjpars,data)

%{ 
    This function solves the time-varying in vivo version of the model to
    steady-state and then calculates 2 steady-state beats. 
    Inputs: 
    adjpars         - vector of adjustable parameters 
    data            - input data structure with data and global parameters
    Outputs: 
    outputs         - structure with all pertinent model outputs to plot 
    rout            - residual vector 
    J               - cost functional 
%} 

%% Initialization  

plot_switch = 0; % 0 = off, 1 = on 

% "undo" log from parameters.m
adjpars = exp(adjpars); 

tspan = data.tspan;  
dt    = data.dt; 

T = data.T; 

ODE_TOL = data.gpars.ODE_TOL;

fixpars = data.fixpars;

%% Parameters 

% Unstressed volumes 
V_SA_u = fixpars(11);
V_SV_u = fixpars(12); 
V_PA_u = fixpars(13); 
V_PV_u = fixpars(14); 

% Compliance 
C_SA = adjpars(1); 

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

%% Solve model 

% Use a while loop to allow model to converge to steady-state 
ndone = 0; 
while ndone == 0 

    opts = odeset('Mass',M,'RelTol',ODE_TOL,'AbsTol',ODE_TOL);
    sol  = ode15s(@model,[tspan(1) tspan(end)],init,opts,adjpars,data);
    
    % Plot intermediate states if plot_switch is "on"
    if plot_switch == 1 
        % Displacements states 1-4
        figure(111)
        clf
        plot(sol.x,sol.y(1:4,:))
        legend('1','2','3','4')
        
        % Sarcomere lengths states 5-7
        figure(112)
        clf
        plot(sol.x,sol.y(5:7,:))
        legend('5','6','7')
        
        % Volumes states 8-11
        figure(113)
        clf
        plot(sol.x,sol.y(8:9,:))
        legend('8','9')
        
        % Volumes states 12-15
        figure(114)
        clf
        plot(sol.x,sol.y(10:13,:))
        legend('10','11','12','13')
    end 

    if sol.x(end) ~= tspan(end) 
        % Check to see if the model solved to then of tspan. If not, set the
        % initial conditions for the next loop at the start of the previous
        % beat and solve for 10 beats 
        t = sol.x(1):dt:sol.x(end); 
        beats = mod(t,T); 
        x = find(round(beats,3) == 0);
        y = find(t(x) <= t(end)); 
        tspan = tspan(1):dt:tspan(x(y(end)));
    
        sols = deval(sol,tspan);
        init  = sols(:,x(y(end-1))); 
    
        tspan = tspan(x(y(end-1))):dt:tspan(x(y(end-1))) + 10*T; 
         
    else 
        % If the model has successfully solved at least 10 beats, then we can
        % assess whether the model has reached steady state 
        
        % Extract sarcomere lengths and systemic arterial pressure (Psa)     
        sols = deval(sol,tspan);
            
        Lsc_LV  = sols(5,:) * 1e6;
        Lsc_SEP = sols(6,:) * 1e6; 
        Lsc_RV  = sols(7,:) * 1e6; 
        P_SA    = sols(9,:) / C_SA * 7.5; 
        
        % Plot extracted quantities if plot_switch is "on" 
        if plot_switch == 1
            figure(101)
            clf
            hold on 
            plot(tspan,Lsc_LV) 
            
            figure(102)
            clf
            hold on 
            plot(tspan,Lsc_SEP) 
            
            figure(103)
            clf
            hold on 
            plot(tspan,Lsc_RV) 
            
            figure(104)
            clf % new fig with each time series
            hold on 
            plot(tspan,P_SA) 
        end
        
        % Find the last 5 beats of the simulation 
        xx = find(tspan >= tspan(end) - 5*T); 
        
        % Set a peak threshold as half of the amplitude of the last 5 beats 
        threshold_LV  = (max(Lsc_LV(xx))  - min(Lsc_LV(xx)))/2;
        threshold_SEP = (max(Lsc_SEP(xx)) - min(Lsc_SEP(xx)))/2;
        threshold_RV  = (max(Lsc_RV(xx))  - min(Lsc_RV(xx)))/2;
        
        % Determine the length of half of the cardiac cycle 
        half_per = round((T/2)/dt); 
        
        % Find peaks for the last 5 beats 
        [pks_Lsc_LV, loc_pks_Lsc_LV]  = findpeaks(...
            Lsc_LV,'MinPeakDistance',half_per,'MinPeakProminence',threshold_LV); 
        [pks_Lsc_SEP,loc_pks_Lsc_SEP] = findpeaks(...
            Lsc_SEP,'MinPeakDistance',half_per,'MinPeakProminence',threshold_SEP); 
        [pks_Lsc_RV, loc_pks_Lsc_RV]  = findpeaks(...
            Lsc_RV,'MinPeakDistance',half_per,'MinPeakProminence',threshold_RV); 
        [pks_P_SA,   loc_pks_P_SA]    = findpeaks(...
            P_SA,'MinPeakDistance',half_per); 
        
        % Exclude the last peak (so there are 4 peaks)
        pks_Lsc_LV  = pks_Lsc_LV(end-5:end-1); 
        pks_Lsc_SEP = pks_Lsc_SEP(end-5:end-1); 
        pks_Lsc_RV  = pks_Lsc_RV(end-5:end-1);
        pks_P_SA    = pks_P_SA(end-5:end-1); 
        
        % Find the locations of the peaks 
        loc_pks_Lsc_LV  = loc_pks_Lsc_LV(end-5:end-1); 
        loc_pks_Lsc_SEP = loc_pks_Lsc_SEP(end-5:end-1); 
        loc_pks_Lsc_RV  = loc_pks_Lsc_RV(end-5:end-1); 
        loc_pks_P_SA    = loc_pks_P_SA(end-5:end-1); 
        
        % Find the times where the peaks occur 
        t_pks_Lsc_LV  = tspan(loc_pks_Lsc_LV);
        t_pks_Lsc_SEP = tspan(loc_pks_Lsc_SEP);
        t_pks_Lsc_RV  = tspan(loc_pks_Lsc_RV);
        t_pks_P_SA    = tspan(loc_pks_P_SA); 
        
        % Create a linear regression through the peaks 
        pf_Lsc_LV  = polyfit(t_pks_Lsc_LV,pks_Lsc_LV,1); 
        pf_Lsc_SEP = polyfit(t_pks_Lsc_SEP,pks_Lsc_SEP,1); 
        pf_Lsc_RV  = polyfit(t_pks_Lsc_RV,pks_Lsc_RV,1); 
        pf_P_SA    = polyfit(t_pks_P_SA,pks_P_SA,1); 
        
        % Extract the slope of the regression line 
        slope_Lsc_LV  = pf_Lsc_LV(1);
        slope_Lsc_SEP = pf_Lsc_SEP(1);
        slope_Lsc_RV  = pf_Lsc_RV(1);
        slope_P_SA    = pf_P_SA(1);
        
        % Plot regression line through peaks if plot_switch is "on" 
        if plot_switch == 1
            % Draw the regression line through the peaks
            y_Lsc_LV  = polyval(pf_Lsc_LV,tspan); 
            y_Lsc_SEP = polyval(pf_Lsc_SEP,tspan); 
            y_Lsc_RV  = polyval(pf_Lsc_RV,tspan); 
            y_P_SA    = polyval(pf_P_SA,tspan); 
            
            % LV sarcomere length 
            figure(101)
            hold on 
            plot(t_pks_Lsc_LV,pks_Lsc_LV,'r*')
            plot(tspan,y_Lsc_LV,'k')
            
            % SEP sarcomere length 
            figure(102)
            hold on 
            plot(t_pks_Lsc_SEP,pks_Lsc_SEP,'r*')
            plot(tspan,y_Lsc_SEP,'k')
            
            % RV sarcomere length 
            figure(103)
            hold on 
            plot(t_pks_Lsc_RV,pks_Lsc_RV,'r*')
            plot(tspan,y_Lsc_RV,'k')
            
            % Systemic arterial pressure 
            figure(104)
            hold on 
            plot(t_pks_P_SA,pks_P_SA,'r*')
            plot(tspan,y_P_SA,'k') 
        end
        
        % If the slope is sufficiently small (i.e. flat), we have reached 
        % steady state 
        slope_lim = 1e-3;
            % Stopping criteria 
            if abs(slope_P_SA) < slope_lim && abs(slope_Lsc_LV) < slope_lim && ...
                    abs(slope_Lsc_SEP) < slope_lim && abs(slope_Lsc_RV) < slope_lim
                ndone = 1; 
            end 
            
        % If we have not reached steady-state, solve the model for 10 more
        % beats and reassess convergence 
        beats = mod(tspan,T); 
        x = find(round(beats,3) == 0);
        y = find(tspan(x) <= tspan(end));
        tspan = tspan(x(y(end-1))):dt:tspan(x(y(end-1))) + 10*T;
        init  = sols(:,x(y(end-1))); 
    end 
end

% After determining that the model is in steady-state, solve 2 more beats 
time = [0:dt:2*T]; 
sol  = ode15s(@model,[time(1) time(end)],init,opts,adjpars,data);
sols = deval(sol,time);
sols = sols'; 

%% Calculate other time-varying model quantities (pressures, flows, etc.) 

o = zeros(38,length(time));  
for i = 1:length(time) 
    [~,o(:,i)] = model(time(i),sols(i,:),adjpars,data);
end 

%% Create output structure  

outputs.time = time; 

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
volumes.V_SA = (sols(:,9) + V_SA_u) * 1e6; 
volumes.V_SV = (sols(:,10) + V_SV_u) * 1e6; 
volumes.V_RV = sols(:,11) * 1e6; 
volumes.V_PA = (sols(:,12) + V_PA_u) * 1e6; 
volumes.V_PV = (sols(:,13) + V_PV_u) * 1e6; 
volumes.Vtot = (sum(sols(end,8:13)) + V_SA_u + V_SV_u + V_PA_u + V_PV_u) * 1e6; 

% Convert kPa to mmHg
pressures.P_LV = o(1,:) * 7.5; 
pressures.P_SA = o(2,:) * 7.5; 
pressures.P_SV = o(3,:) * 7.5; 
pressures.P_RV = o(4,:) * 7.5; 
pressures.P_PA = o(5,:) * 7.5; 
pressures.P_PV = o(6,:) * 7.5; 

% Convert m^3 to cm^3
wallvolumes.Vm_LV  = o(7,:)  * 1e6; 
wallvolumes.Vm_SEP = o(8,:) * 1e6; 
wallvolumes.Vm_RV  = o(9,:) * 1e6; 

% Convert m^2 to cm^2
areas.Am_LV  = o(10,:) * 1e4; 
areas.Am_SEP = o(11,:) * 1e4; 
areas.Am_RV  = o(12,:) * 1e4; 

% Convert m^(-1) to cm^(-1)
curvatures.Cm_LV  = o(13,:) * 1e-2;
curvatures.Cm_SEP = o(14,:) * 1e-2;
curvatures.Cm_RV  = o(15,:) * 1e-2; 

strains.eps_LV  = o(16,:); 
strains.eps_SEP = o(17,:); 
strains.eps_RV  = o(18,:); 

stresses.passive.sigma_pas_LV  = o(19,:);
stresses.passive.sigma_pas_SEP = o(20,:);
stresses.passive.sigma_pas_RV  = o(21,:);

stresses.active.sigma_act_LV  = o(22,:);
stresses.active.sigma_act_SEP = o(23,:);
stresses.active.sigma_act_RV  = o(24,:);

stresses.total.sigma_LV  = o(25,:);
stresses.total.sigma_SEP = o(26,:);
stresses.total.sigma_RV  = o(27,:);

% Convert m^3 s^(-1) to L min^(-1)
flows.Q_m_valve = o(28,:) * 1e3 * 60; 
flows.Q_a_valve = o(29,:) * 1e3 * 60; 
flows.Q_t_valve = o(30,:) * 1e3 * 60; 
flows.Q_p_valve = o(31,:) * 1e3 * 60; 

flows.Q_SA = o(32,:) * 1e3 * 60; 
flows.Q_PA = o(33,:) * 1e3 * 60; 

tensions.Tm_LV = o(34,:);
tensions.Tm_SEP = o(35,:); 
tensions.Tm_RV = o(36,:); 

activation.y_v = o(37,:)';

pressures.P_peri = o(38,:) * 7.5'; 

% Q_m_valve = flows.Q_m_valve; 
% Q_a_valve = flows.Q_a_valve; 
% 
% i_ES = find(diff(Q_m_valve) > 0,1,'first'); 
% i_ED = find(diff(Q_a_valve) > 0,1,'first'); 
%     
% ES = time(i_ES); 
% ED = time(i_ED); 
% 
% timepoints.ES = ES; 
% timepoints.ED = ED; 
% outputs.timepoints    = timepoints; 

outputs.volumes       = volumes; 
outputs.pressures     = pressures; 
outputs.displacements = displacements; 
outputs.areas         = areas;
outputs.wallvolumes   = wallvolumes; 
outputs.curvatures    = curvatures; 
outputs.strains       = strains; 
outputs.stresses      = stresses;
outputs.lengths       = lengths; 
outputs.flows         = flows; 
outputs.tensions      = tensions; 
outputs.activation    = activation; 

%% Sensitivity Analysis 

rout = []; 

% catch 
%        outputs = [];
%        rout = 100;
% end

J = rout'*rout;

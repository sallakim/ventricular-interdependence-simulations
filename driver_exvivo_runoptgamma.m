% This script optimizes gamma

clear all
%close all
tic

addpath ExVivo
addpath Core
addpath OPT_codes

warning('on')
                                
%% Load pseudodata

data = makedatastructure;

a_eta_Vtot = [.6:.1:2];
data.a_eta_Vtot = a_eta_Vtot; 

%% Get nominal parameter values

[pars,UB,LB,data] = parameters(data); 

INDMAP = 23; % gamma

ALLPARS  = pars;
ODE_TOL  = 1e-8; 
DIFF_INC = sqrt(ODE_TOL);

gpars.INDMAP   = INDMAP;
gpars.ALLPARS  = ALLPARS;
gpars.ODE_TOL  = ODE_TOL;
gpars.DIFF_INC = DIFF_INC;

data.gpars = gpars;

%% Optimization - fmincon

optx   = pars(INDMAP); 
opthi  = UB(INDMAP);
optlow = LB(INDMAP);

maxiter = 40; 
mode    = 2; 
nu0     = 2.d-1; 

[xopt, histout, costdata, jachist, xhist, rout, sc] = ...
     newlsq_v2(optx,'opt_wrap',1.d-4,maxiter,...
     mode,nu0,opthi,optlow,data); 

pars_opt = pars;
pars_opt(INDMAP) = xopt; 

% run model with optimized parameters 
data.gamma_opt = exp(xopt); 
[adjpars,~,~,data] = parameters(data); 
[outputs,rout,J] = model_sol_gammaopt(adjpars,data);

optpars = exp(pars_opt);
disp('optimized gamma')
disp([INDMAP' optpars(INDMAP)])

save opt_gamma.mat 

elapsed_time = toc;
elapsed_time = elapsed_time/60

%% Plot figures 

a_EDP_LV = outputs.EDPVRs.a_EDP_LV; 
a_EDV_LV = outputs.EDPVRs.a_EDV_LV; 
a_EDP_RV = outputs.EDPVRs.a_EDP_RV; 
a_EDV_RV = outputs.EDPVRs.a_EDV_RV; 

V_LV_EDPVR  = outputs.EDPVRs.V_LV_EDPVR;
V_RV_EDPVR  = outputs.EDPVRs.V_RV_EDPVR;
P_LV_EDPVR  = outputs.EDPVRs.P_LV_EDPVR;
P_RV_EDPVR  = outputs.EDPVRs.P_RV_EDPVR;

hfig201 = figure(201); 
clf
hold on 
h1 = plot(V_LV_EDPVR, P_LV_EDPVR,'r','linewidth',2); 
h2 = plot(a_EDV_LV, a_EDP_LV, 'ro', 'markersize',10); 
legend([h1 h2],'Klotz','Model','location','northwest')
set(gca,'FontSize',20)
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
ylim([0 60])


hfig202 = figure(202); 
clf
hold on 
h1 = plot(V_RV_EDPVR, P_RV_EDPVR,'b','linewidth',2); 
h2 = plot(a_EDV_RV, a_EDP_RV, 'bo', 'markersize',10); 
legend([h1 h2],'Klotz','Model','location','northwest')
set(gca,'FontSize',20)
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
ylim([0 60])


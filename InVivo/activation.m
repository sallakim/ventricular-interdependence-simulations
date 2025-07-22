function y = activation(t,data)

%% Unpack data structure

HR = data.HR; 

k_TS = data.fixpars(15);
k_TR = data.fixpars(16); 

%% Compute Activation Function

% Heart period (s) 
T = 60/HR; 

% Time to maximal systole 
TS_v = k_TS * T; 

% Time from maximal systole to relaxation 
TR_v = k_TR * T; 


tc = mod(t,T);
if tc >= 0 && tc < TS_v
    y = 0.5*(1 - cos(pi*(tc)/TS_v)); 
elseif tc >= TS_v && tc < TR_v + TS_v 
    y = 0.5*(1 + cos(pi*(tc - TS_v)/TR_v)); 
else
    y = 0; 
end 

end 

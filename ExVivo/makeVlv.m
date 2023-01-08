


t = 0:.01:1; 

EDV_lv = 125; 
ESV_lv = 50; 

q_down = 25; 
q_up = 50; 

s_down = .25;
s_up = .7;

Vlv_d = (EDV_lv - ESV_lv) * 1./ (1 + exp(q_down * (t - s_down))) + ESV_lv;

Vlv_u = (EDV_lv - ESV_lv) * 1./ (1 + exp(-q_up * (t - s_up))) + ESV_lv; 

Vlv = Vlv_d + Vlv_u - ESV_lv; 

figure(1)
clf
plot(t,Vlv_d,t,Vlv_u,t,Vlv)
% visualize the time course
clear;
xy = load('TriSegAnim_normal.mat').xy;
% xy = [t, xm_lv, xm_sep,  xm_rv, ym, t_lv, t_sep,  t_rv]; % data structure
% format
slides = [1, 500, 1000];
t = xy(slides , 1);
xm_lv = xy(slides , 2);


ym = xy(slides ,5);

% Circle Function For Angles In Degrees
plotArc = @(radius,deg_ang)  [radius*cosd(deg_ang);  radius*sind(deg_ang)]; 

Vw_lv = 8.2;
% Vw_sep = 4.1;
% Vw_rv = 2.6;

Cm=2*xm_lv./max(xm_lv.^2 + ym.^2, 1e-6);
R = 1./Cm;

x0 = xm_lv - R;
a = atan2(ym, x0);

a_start = -pi/2 -a;
a_pos = xm_lv > 0;
a_start(a_pos)= -pi/2 + a(a_pos);

a_stop = -pi/2 + a;
a_stop(xm_lv > 0) =  3*pi/2 - a(a_pos);

% Cm_sep=2*xm_sep/max(xm_sep^2 + ym^2, 1e-6);
% Cm_rv=2*xm_rv/max(xm_rv^2 + ym^2, 1e-6);
plotArc = @(radius,arc)  [radius.*cos(arc);  radius.*sin(arc)]; 

i = 1;
arc = linspace(a_start(i), a_stop(i), 20);
xy = plotArc(R(i), arc);
plot(xy(1, :), xy(2, :))



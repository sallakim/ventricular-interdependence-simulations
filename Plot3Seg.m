% visualize the time course
clear;
xy = load('TriSegAnim_normal.mat').xy;
% xy = [t, xm_lv, xm_sep,  xm_rv, ym, t_lv, t_sep,  t_rv]; % data structure
% format
slides = [1:10:900];
t = xy(slides , 1);
xm_lv = xy(slides , 2);
xm_sep = xy(slides , 3);
xm_rv = xy(slides , 4);
ym = xy(slides ,5);

% Circle Function For Angles In Degrees
plotArc = @(radius,deg_ang)  [radius*cosd(deg_ang);  radius*sind(deg_ang)]; 

Vw_lv = 8.2;
Vw_sep = 4.1;
Vw_rv = 2.6;

Cm_lv=2*xm_lv./max(xm_lv.^2 + ym.^2, 1e-6);
Cm_sep=2*xm_sep./max(xm_sep.^2 + ym.^2, 1e-6);
Cm_rv=2*xm_rv./max(xm_rv.^2 + ym.^2, 1e-6);

R_lv = 1./Cm_lv;
R_sep = 1./Cm_sep;
R_rv = 1./Cm_rv;

x0_lv = xm_lv - R_lv;
x0_sep = xm_sep - R_sep;
x0_rv = xm_rv - R_rv;

a_lv = atan2(ym, x0_lv);
a_sep = atan2(ym, -x0_sep);
a_rv = atan2(ym, -x0_rv);

a_start_lv = 0*pi/2 - a_lv;
a_stop_lv = 0*pi/2 + a_lv;

a_start_sep = 0*pi/2 -a_sep;
a_stop_sep = 0*pi/2 + a_sep;

a_start_rv = 0*pi/2 -a_rv;
a_stop_rv = 0*pi/2 + a_rv;

plotArc = @(x0, radius,arc)  [radius.*cos(arc) + x0;  radius.*sin(arc)]; 

i = 1;
figure(101);clf; hold on;
for i = 1:length(t)

    arc_lv = linspace(a_start_lv(i), a_stop_lv(i), 20);
    xy_lv = plotArc(x0_lv(i), R_lv(i), arc_lv);
    plot(xy_lv(1, :), xy_lv(2, :))
    
    arc_sep = linspace(a_start_sep(i), a_stop_sep(i), 20);
    xy_sep = plotArc(x0_sep(i), R_sep(i), arc_sep);
    plot(xy_sep(1, :), xy_sep(2, :))
    
    arc_rv = linspace(a_start_rv(i), a_stop_rv(i), 20);
    xy_rv = plotArc(x0_rv(i), R_rv(i), arc_rv);
    plot(xy_rv(1, :), xy_rv(2, :))

end
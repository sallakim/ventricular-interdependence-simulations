% visualize the time course
clear;
figure(102);clf; hold on;

defpos = [200          69        1199         922];

files = {'TriSegAnim_normal.mat', ...
         'TriSegAnim_LVSD.mat', ...
         'TriSegAnim_LVDD.mat', ...
         'TriSegAnim_RVSD.mat', ...
         'TriSegAnim_RVDD.mat'};
titles = {'Healthy normal', ...
          'Left ventricular Systolic dysfunction', ...
          'Left ventricular Diastolic dysfunction', ...
          'Right ventricular Systolic dysfunction',...
          'Right ventricular Diastolic dysfunction'};

for caseIter = 1:length(files)
    xy = load(['..\' files{caseIter}]).xy;
    
    set(gcf, 'Position', defpos);
    % xy = [t, xm_lv, xm_sep,  xm_rv, ym, t_lv, t_sep,  t_rv]; % data structure
    % format
    % slides = round(linspace(1, 800, 8)/10)*10;
    slides = [97, 374, 556, 1000];
    slides_titles = {"Start systole", "End systole", "Start diastole", "End diastole"};
    t = xy(slides, 1)';
    xm_lv = xy(slides , 2)';
    xm_sep = xy(slides , 3)';
    xm_rv = xy(slides , 4)';
    ym = xy(slides ,5)';
    
    Vw_lv = 82; % cm3
    Vw_sep = 41; % cm3
    Vw_rv = 26; % cm3
    
    % Calculation the wall thickness
    VmInner = @(h, xm, ym) 1/6*pi*(xm + h).*(3*ym.^2 + (xm + h).^2);
    VmOuter = @(h, xm, ym) 1/6*pi*(xm - h).*(3*ym.^2 + (xm - h).^2);
    % this should be zero, as Vw = VmOUter - VmInner;
    VwBalance = @(h, xm, ym, Vw) -VmOuter(h, xm, ym) + VmInner(h, xm, ym) - Vw; 
    for i = 1:length(xm_lv)
        h_lv(i) = fzero(@(h)VwBalance(h, xm_lv(i), ym(i), Vw_lv),  [-10 10]);
        h_sep(i) = fzero(@(h)VwBalance(h, xm_sep(i), ym(i), Vw_sep),  [-10 10]);
        h_rv(i) = fzero(@(h)VwBalance(h, xm_rv(i), ym(i), Vw_rv),  [-10 10]);
    end
    % %add up to the thicknesses
    % % inner lining
    % xm_lv = xm_lv + h_lv;
    % xm_sep = xm_sep - h_sep;
    % xm_rv = xm_rv - h_rv;
    % % ym = ym - h_lv;
    % ym_lv = ym - h_lv;
    % ym_sep = ym;
    % ym_rv = ym - h_rv;
    
    % % outer lining
    % xm_lv = xm_lv - h_lv;
    % xm_sep = xm_sep + h_sep;
    % xm_rv = xm_rv + h_rv;
    % ym_lv = ym + h_lv;
    % ym_sep = ym;
    % ym_rv = ym + h_rv;
    
    
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
    co = colororder;
    for i = 1:length(slides)
        subplot(length(files), length(slides), i + length(slides)*(caseIter-1));hold on;
        set(gca, "FontSize", 14);
axis([-5 5 -5.5 5.5])
axis equal;

        if i == 1
            % only first have ticks and label
            ylabel("cm");
            set(gca, "YTick", [-5, 0, 5]);
            set(gca, "YTicklabel", [-5, 0, 5], 'TickLength', [0.0250 0.0250]);
        else
            % no ticks for insides
            set(gca, "YTick", []);
            set(gca, "YTicklabel", []);
        end
        if length(files) == caseIter
            xlabel({"cm", slides_titles{i}});
            set(gca, "XTick", [-5, 0, 5]);
            set(gca, "XTicklabel", [-5, 0, 5], 'TickLength', [0.0250 0.0250]);
        else
% no ticks for insides
            set(gca, "XTick", []);
            set(gca, "XTicklabel", []);
        end
        if i == 1% round(length(slides)/2)
            title(titles{caseIter}, "HorizontalAlignment","left");
        end
        c = co(i, :);
        % Outer LV
        arc_lv = linspace(a_start_lv(i), a_stop_lv(i), 20);
        xy_lv = plotArc(x0_lv(i), R_lv(i) - h_lv(i), arc_lv);
        plot(xy_lv(1, :), xy_lv(2, :), 'LineWidth',2, 'Color', c);
        % midline LV
        arc_lv = linspace(a_start_lv(i), a_stop_lv(i), 20);
        xy_lv = plotArc(x0_lv(i), R_lv(i), arc_lv);
        plot(xy_lv(1, :), xy_lv(2, :), '--', 'LineWidth',0.5, 'Color', c);    
        % inner LV
        arc_lv = linspace(a_start_lv(i), a_stop_lv(i), 20);
        xy_lv = plotArc(x0_lv(i), R_lv(i) + h_lv(i), arc_lv);
        plot(xy_lv(1, :), xy_lv(2, :), 'LineWidth',2, 'Color', c)
    
    
        % LV side
        arc_sep = linspace(a_start_sep(i), a_stop_sep(i), 20);
        xy_sep = plotArc(x0_sep(i), R_sep(i) - h_sep(i), arc_sep);
        plot(xy_sep(1, :), xy_sep(2, :), 'LineWidth',2, 'Color', c)
        % mid SEP
        arc_sep = linspace(a_start_sep(i), a_stop_sep(i), 20);
        xy_sep = plotArc(x0_sep(i), R_sep(i), arc_sep);
        plot(xy_sep(1, :), xy_sep(2, :), '--', 'Color', c)
        % RV side SEP
        arc_sep = linspace(a_start_sep(i), a_stop_sep(i), 20);
        xy_sep = plotArc(x0_sep(i), R_sep(i) + h_sep(i), arc_sep);
        plot(xy_sep(1, :), xy_sep(2, :), 'LineWidth',2, 'Color', c)
        
        % RV inner
        arc_rv = linspace(a_start_rv(i), a_stop_rv(i), 20);
        xy_rv = plotArc(x0_rv(i), R_rv(i) - h_rv(i), arc_rv);
        plot(xy_rv(1, :), xy_rv(2, :), 'LineWidth',2, 'Color', c)
        % mid
        arc_rv = linspace(a_start_rv(i), a_stop_rv(i), 20);
        xy_rv = plotArc(x0_rv(i), R_rv(i), arc_rv);
        plot(xy_rv(1, :), xy_rv(2, :), '--', 'Color', c)
        % RV outer
        arc_rv = linspace(a_start_rv(i), a_stop_rv(i), 20);
        xy_rv = plotArc(x0_rv(i), R_rv(i) + h_rv(i), arc_rv);
        plot(xy_rv(1, :), xy_rv(2, :), 'LineWidth',2, 'Color', c)
    
        % xlim([-5 5]);ylim([-5 5]);axis equal;
        % axis equal;
    end
end 

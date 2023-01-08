function TriSegAnimExport(sim_outputs, exp_name)
  
% xy0 = [0, -4.6627, 2.90348, 6.26344,3.50013, 1, 0.5, 0];
    
%   dt = 0.1;
%   endtime = 10;
%   
%   xy = xy0;
%   for t = dt:dt:endtime
%       
%      
%      xy = [xy;[t, xy0(2) + sin(t), xy0(3) + cos(t), xy0(4) + sin(t),xy0(5) + t/5], xy0(6),xy0(7),xy0(8)];
%   end
  
% %   not needed, just shos the structure
%   xyt = array2table({'time', 'xm_LV', 'xm_SEP', 'xm_RV', 'ym', 'T LV', 'T SEP', 'T RV'});
%   plotFigs = false;
%   
%   if plotFigs
%       figure(1);clf;
%         subplot(211); title("Positions")
%         plot(xy(:, 1), xy(:, 2),xy(:, 1), xy(:, 3), xy(:, 1), xy(:, 4),xy(:, 1), xy(:, 5), 'LineWidth', 2);
%         legend('[2] xm_{LV}', '[3] xm_{SEP}', '[4] xm_{RV}', '[5] ym')
%         subplot(212); title("colouring (Tensions), constant by default");
%         plot(xy(:, 1), xy(:, 6),xy(:, 1), xy(:, 7), xy(:, 1), xy(:, 8), 'LineWidth', 3);
%       legend('[6] T LV', '[7] T SEP', '[4] T RV');
%   end

xy = [sim_outputs.time', ...
    -sim_outputs.displacements.xm_LV,...
    sim_outputs.displacements.xm_SEP,...
    sim_outputs.displacements.xm_RV,...
    sim_outputs.displacements.ym,...
    sim_outputs.tensions.Tm_LV',...
    sim_outputs.tensions.Tm_SEP',...
    sim_outputs.tensions.Tm_RV'...
    ];
    
  
    
    save(['TriSegAnim_' exp_name '.mat'], 'xy');
    
    
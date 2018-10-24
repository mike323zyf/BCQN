function [ ] = plot_result( u_x, frame )
%PLOT_RESULT Summary of this function goes here
%   Detailed explanation goes here

global obj_tri ver_num cmin cmax K q

q_x = K * u_x + q;

color = energy_color(q_x);

r_q_x = reshape(q_x, 2, ver_num)';

TR = triangulation(obj_tri, r_q_x);
indB = TR.freeBoundary();
intB = [indB(:, 1); indB(1, 1)];
boundary_line_width = 1;

file = figure(1);

clf

minmin = min(r_q_x);
maxmax = max(r_q_x);

min_x = minmin(1);
min_y = minmin(2);
max_x = maxmax(1);
max_y = maxmax(2);

cen_x = 0.5 * (min_x + max_x);
cen_y = 0.5 * (min_y + max_y);

radius = max(max_x - min_x, max_y - min_y);
radius = radius * 0.6;

axis([cen_x - radius cen_x + radius cen_y - radius cen_y + radius]);
caxis([cmin cmin + cmax])



patch('Faces', obj_tri, 'Vertices', r_q_x, 'FaceVertexCData', color(:, 1), 'FaceColor','flat', 'EdgeAlpha',0.1);
% patch('Faces', obj_tri, 'Vertices', r_q_x, 'FaceVertexCData', color(:, 1), 'FaceColor','flat', 'EdgeAlpha',0.0);

% colorbar;
colormap(getDistortionColormap())

line(r_q_x(intB, 1, end), r_q_x(intB, 2, end), 'color', 'b', 'linewidth', boundary_line_width);

axis square
% axis off
% 
title(frame);
% 
drawnow;
% 
% getframe;
% % 
% set(gca,'units','centimeters')
% set(gca,'Position',[1.8 0.8 14.8 14.8]);
% set(gcf, 'PaperUnits','centimeters');
% set(gcf, 'PaperSize', [12 12]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition',[0 0 12 12]);
% 
% print(file, '-dbmp', sprintf(strcat('../result/image%d.bmp'), frame), '-r400')  

end


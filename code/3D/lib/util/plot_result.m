function [ ] = plot_result( u_x, frame )
%PLOT_RESULT Summary of this function goes here
%   Detailed explanation goes here

global tet_mesh ver_num cmin cmax K q

obj_tet = tet_mesh.t;

q_x = K * u_x + q;

color = energy_color(q_x);

r_q_x = reshape(q_x, 3, ver_num)';

file = figure(1);

% set(gcf, 'Render', 'opengl');

clf

minmin = min(r_q_x);
maxmax = max(r_q_x);

min_x = minmin(1);
min_y = minmin(2);
min_z = minmin(3);
max_x = maxmax(1);
max_y = maxmax(2);
max_z = maxmax(3);

cen_x = 0.5 * (min_x + max_x);
cen_y = 0.5 * (min_y + max_y);
cen_z = 0.5 * (min_z + max_z);

radius = max(max(max_x - min_x, max_y - min_y), max_z - min_z);
radius = max(radius, 0.1) * 0.5 * 1.1;

tetramesh(obj_tet, r_q_x);
% patch('Faces', obj_tri, 'Vertices', r_q_x, 'FaceVertexCData', color(:, 1), 'FaceColor','flat');

axis([cen_x - radius cen_x + radius cen_y - radius cen_y + radius cen_z - radius cen_z + radius]);
% caxis([cmin cmin + cmax])

axis square

colorbar;

title(frame);

drawnow;


% getframe;
% 
% set(gca,'units','centimeters')
% set(gca,'Position',[0.8 0.8 14.8 14.8]);
% set(gcf, 'PaperUnits','centimeters');
% set(gcf, 'PaperSize', [12 12]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition',[0 0 12 12]);
% 
% print(file, '-dbmp', sprintf(strcat('../result/image%d.bmp'), frame), '-r300') 


end


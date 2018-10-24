
global c1 c2 d1 tol_x_cnt tol_f_cnt stop_cnt tol_x tol_f amips_s USE_TBB

% c1 = 12;
% c2 = 4;
% d1 = 18;
% Neo-Hookean (~Mooney-Rivlin for 2D)

scale = 1;

c1 = 8 * scale; % shear
d1 = 2 * scale;  % bulk

c2 = 0; % meaningless in 2D 

tol_x_cnt = 0;
tol_f_cnt = 0;

tol_x = 1e-10;
tol_f = 1e-6;

stop_cnt = 5;

amips_s = 1;

% USE_TBB = 1;
USE_TBB = 0;


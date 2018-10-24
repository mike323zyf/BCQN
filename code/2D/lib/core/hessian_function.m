function [ H ] = hessian_function( u_x, clamp )
%HESSIAN_FUNCTION Summary of this function goes here
%   Detailed explanation goes here

global tri_num X_g_inv tri_areas exam_mesh amips_s F_dot ver_num K q c1 c2 d1

q_x = K * u_x + q;

[rows, cols, vals] = energy_hessian_mex(tri_num, X_g_inv, tri_areas, exam_mesh.f, q_x, get_energy_type(), amips_s, F_dot, ver_num, c1, c2, d1, clamp);

rows = rows + 1;
cols = cols + 1;

H = sparse(rows, cols, vals, 2 * ver_num, 2 * ver_num);

H = K' * H * K;

% H = reshape(T, 2 * ver_num, 2 * ver_num)';
% 
% H = K' * H * K;

end


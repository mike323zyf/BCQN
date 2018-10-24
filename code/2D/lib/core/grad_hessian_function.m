function [ grad, H ] = grad_hessian_function( u_x, clamp )
%GRAD_HESSIAN_FUNCTION Summary of this function goes here
%   Detailed explanation goes here

global tri_num X_g_inv tri_areas exam_mesh amips_s F_dot ver_num K q c1 c2 d1 x2u J_u J_ui JT_u JT_ui JV_u JTV_u wu bu

q_x = K * u_x + q;

JTJ_in = [J_ui(1), J_ui(2), JT_ui(1), JT_ui(2)]';

[J_q, JV_u, JTV_u, wu, bu, rows, cols, vals] = energy_hessian_mex(tri_num, X_g_inv, tri_areas, exam_mesh.f, q_x, get_energy_type(), amips_s, F_dot, ver_num, c1, c2, d1, clamp, J_u, JT_u, JTJ_in, x2u);

grad = K' * J_q;

rows = rows + 1;
cols = cols + 1;

H = sparse(rows, cols, vals, 2 * ver_num, 2 * ver_num);

H = K' * H * K;

end


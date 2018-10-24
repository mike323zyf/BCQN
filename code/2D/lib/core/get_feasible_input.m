function [ output ] = get_feasible_input()
%GET_FEASIBLE_INPUT Summary of this function goes here
%   Detailed explanation goes here

global dirichlet q_rest ver_num exam_mesh tri_num

ver = exam_mesh.v(:, 1:2);
ver = reshape(ver', 2 * ver_num, 1);

[X_g_inv_m, tri_areas_m, F_dot_m, rows, cols, vals, x2u, J_u, J_ui, JT_u, JT_ui] = precompute_mex(tri_num, exam_mesh.f, ver_num, ver, 0, dirichlet);

rows = rows + 1;
cols = cols + 1;

A = sparse(rows, cols, vals, 2 * ver_num, 2 * ver_num);

A = 2 * A;

dir_i = find(dirichlet > 0);

c_num = size(dir_i, 1);

U = sparse(c_num, 2 * ver_num);

for i = 1 : c_num
   
    U(i, dir_i(i)) = 1;
    
end

KKT_A = sparse(2 * ver_num + c_num, 2 * ver_num + c_num);
KKT_A(1 : 2 * ver_num, 2 * ver_num + 1 : 2 * ver_num + c_num) = U';
KKT_A(2 * ver_num + 1 : 2 * ver_num + c_num, 1 : 2 * ver_num) = U;
KKT_b = zeros(2 * ver_num + c_num, 1);

KKT_b(2 * ver_num + 1 : 2 * ver_num + c_num) = q_rest(dir_i);
KKT_A(1 : 2 * ver_num, 1 : 2 * ver_num) = A;


KKT_x = KKT_A \ KKT_b;

output = KKT_x(1 : 2 * ver_num);

end


function [ H ] = get_precomputed_split_info( is_uv_mesh )
%GET_PRECOMPUTED_INFO Summary of this function goes here
%   Detailed explanation goes here

global tri_num ver_num exam_mesh  X_g_inv tri_areas F_dot K total_volume dirichlet x2u J_u J_ui JT_u JT_ui perimeter_norm Ks free_num order

ver = 0;

if is_uv_mesh == 1
    
    ver = exam_mesh.v(:, 1:3);

    ver = reshape(ver', 3 * ver_num, 1);
    
else
    
    ver = exam_mesh.v(:, 1:2);

    ver = reshape(ver', 2 * ver_num, 1);
    
end

[X_g_inv_m, tri_areas_m, F_dot_m, rows, cols, vals, x2u, J_u, J_ui, JT_u, JT_ui, perimeter] = precompute_mex(tri_num, exam_mesh.f, ver_num, ver, is_uv_mesh, dirichlet);

% check = norm(perimeter)
perimeter_norm = norm(perimeter) * get_energy_characteristic()

X_g_inv = reshape(X_g_inv_m, tri_num, 2, 2);
tri_areas = reshape(tri_areas_m, tri_num, 1);
F_dot = reshape(F_dot_m, tri_num, 6, 2, 2);

total_volume = sum(tri_areas);

rows = rows + 1;
cols = cols + 1;

A = sparse(rows, cols, vals, 2 * ver_num, 2 * ver_num);

A = 2 * A;

Ap = K' * A * K;

order = symamd(Ap);

As = Ap(order, order);

free_num = size(K, 2) / 2;

H = As([1:free_num],[1:free_num]);

end


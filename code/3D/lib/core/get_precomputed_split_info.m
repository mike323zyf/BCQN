function [ H ] = get_precomputed_split_info()
%GET_PRECOMPUTED_INFO Summary of this function goes here
%   Detailed explanation goes here

global tet_num ver_num tet_mesh  X_g_inv tet_vols F_dot K total_volume dirichlet x2u J_u J_ui JT_u JT_ui perimeter_norm free_num order


ver = tet_mesh.v(:, :);
ver = reshape(ver', 3 * ver_num, 1);

[X_g_inv_m, tet_vols_m, F_dot_m, rows, cols, vals, x2u, J_u, J_ui, JT_u, JT_ui, perimeter] = precompute_mex(tet_num, tet_mesh.t, ver_num, ver, dirichlet);

perimeter_norm = norm(perimeter) * get_energy_characteristic();

X_g_inv = reshape(X_g_inv_m, tet_num, 3, 3);
tet_vols = reshape(tet_vols_m, tet_num, 1);
F_dot = reshape(F_dot_m, tet_num, 12, 3, 3);

total_volume = sum(tet_vols);

rows = rows + 1;
cols = cols + 1;

A = sparse(rows, cols, vals, 3 * ver_num, 3 * ver_num);

A = 2 * A;

Ap = K' * A * K;

order = symamd(Ap);

As = Ap(order, order);

free_num = size(K, 2) / 3;

H = As([1:free_num],[1:free_num]);

end


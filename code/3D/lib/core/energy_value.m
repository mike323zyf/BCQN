function [ output ] = energy_value( u )
%ENERGY_VALUE Summary of this function goes here
%   Detailed explanation goes here

global tet_num X_g_inv tet_vols tet_mesh amips_s K q c1 c2 d1

q_target = K * u + q;

output = energy_value_mex(tet_num, X_g_inv, tet_vols, tet_mesh.t, q_target, get_energy_type(), amips_s, c1, c2, d1);

end


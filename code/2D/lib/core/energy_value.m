function [ output ] = energy_value( u )
%ENERGY_VALUE Summary of this function goes here
%   Detailed explanation goes here

global tri_num X_g_inv tri_areas exam_mesh amips_s K q c1 c2 d1

q_target = K * u + q;

output = energy_value_mex(tri_num, X_g_inv, tri_areas, exam_mesh.f, q_target, get_energy_type(), amips_s, c1, c2, d1);

end


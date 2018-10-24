function [ color ] = energy_color( q_target )
%ENERGY_COLOR Summary of this function goes here
%   Detailed explanation goes here

global tri_num X_g_inv tri_areas exam_mesh amips_s c1 c2 d1

color = energy_color_mex(tri_num, X_g_inv, tri_areas, exam_mesh.f, q_target, get_energy_type(), amips_s, c1, c2, d1);

end


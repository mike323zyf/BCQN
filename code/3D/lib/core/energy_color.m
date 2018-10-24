function [ color ] = energy_color( q_target )
%ENERGY_COLOR Summary of this function goes here
%   Detailed explanation goes here

global tet_num X_g_inv tet_vols tet_mesh amips_s c1 c2 d1

color = energy_color_mex(tet_num, X_g_inv, tet_vols, tet_mesh.t, q_target, get_energy_type(), amips_s, c1, c2, d1);

end


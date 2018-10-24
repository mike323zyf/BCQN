function [ ] = save_mesh( u, ver_dim, tri_dim, frame )
%SAVE_MESH Summary of this function goes here
%   Detailed explanation goes here

global K q obj_tri ver_num tri_num

q_x = K * u + q;

save_mesh_mex(tri_num, ver_num, obj_tri, q_x, tri_dim, ver_dim, strcat('../result/frame', int2str(frame), '.obj'));

end


function [ ] = save_mesh( u, frame )
%SAVE_MESH Summary of this function goes here
%   Detailed explanation goes here

global K q tet_mesh ver_num tet_num

q_x = K * u + q;

save_tet_mex(tet_num, ver_num, tet_mesh.t, q_x, strcat('../result/frame', int2str(frame), '.tet'));

end


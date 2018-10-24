function [] = load_input_mesh( file_name )
%LOAD_INPUT_MESH Summary of this function goes here
%   Detailed explanation goes here

global tet_mesh ver_num tet_num

addpath('./data');

tet_mesh = struct('v', [], 't', []);

[vers, tets] = load_tet_mex(strcat(file_name, '.tet'));

tets = tets + 1;

ver_num = size(vers, 1) / 3;
tet_num = size(tets, 1) / 4;

tet_mesh.v = reshape(vers, ver_num, 3);
tet_mesh.t = reshape(tets, tet_num, 4);

end


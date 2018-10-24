function [] = load_input_mesh( file_name )
%LOAD_INPUT_MESH Summary of this function goes here
%   Detailed explanation goes here

global exam_mesh ver_num tri_num obj_tri

addpath('./data');

exam_mesh = struct('v', [], 'n', [], 'u', [], 'f', [], 'e', [], 'bounds', []);

[v, f, uv] = load_mesh_mex(strcat(file_name, '.obj'));

f = f + 1;

ver_num = size(v, 1) / 3;
tri_num = size(f, 1) / 3;
uv_num = size(uv, 1) / 2;

exam_mesh.v = reshape(v, ver_num, 3);
exam_mesh.f = reshape(f, tri_num, 3);

if uv_num > 0
    
    exam_mesh.u = reshape(uv, uv_num, 2);
    
end

obj_tri = exam_mesh.f;

end


function [ V0 ] = get_uv_init( )
%GET_UV_INIT Summary of this function goes here
%   Detailed explanation goes here

global exam_mesh tri_num X_g_inv ver_num

[V0, inds_bF] = computeTutte(exam_mesh.f,exam_mesh.v,true);

q0 = reshape(V0', ver_num * 2, 1);

check = local_injectivity_check_mex(tri_num, exam_mesh.f, q0);

if check == 0
    
    V0 = tutte_embedding_mex(exam_mesh.f, exam_mesh.v);
    
    'Panozo'
    
end

V0 = tutte_embedding_mex(exam_mesh.f, exam_mesh.v);

end


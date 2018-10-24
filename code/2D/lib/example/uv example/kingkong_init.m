
global ver_num tri_num obj_tri exam_mesh dirichlet en_type tol_x_cnt tol_f_cnt stop_cnt tol_x tol_f amips_s cmin cmax K q

colormap(jet);

en_type = 'iso';

load '../data/kingkong.mat';

exam_mesh = struct('v', [], 'n', [], 'u', [], 'f', [], 'e', [], 'bounds', []);

ver_num = size(mesh.v, 1);
tri_num = size(mesh.f, 1);

exam_mesh.v = mesh.v;
exam_mesh.f = mesh.f;
exam_mesh.u = mesh.u;

obj_tri = exam_mesh.f;

% for i = 1 : exam_mesh_num
%    
%     exam_mesh(i)
%     
% end


% precompute_step;

% exam_mesh(1).v(:, 3) = exam_mesh(1).v(:, 3) * 0;
% exam_mesh(1).v(:, 1:2) = exam_mesh(1).u(:, :);
% 
% writeMesh(exam_mesh(1), 'test.obj')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q = reshape(exam_mesh.u(:, 1:2)', ver_num * 2, 1);



mask = [1];

% mask = zeros(ver_num, 1);

% dirichlet = zeros(2 * ver_num, 1);
% 
% for i = 1 : ver_num
%     
%     bd_v = find(exam_mesh.e(i, :) == 1);
%     
%     if size(bd_v, 2) > 0
%        
% %         boundary_ver(i) = 1;
%         dirichlet(2 * i - 1) = 1;
%         dirichlet(2 * i) = 1;
%         
%     end
%     
% end
% 
% q = get_feasible_input();

dirichlet = zeros(2 * ver_num, 1);
dirichlet([2 * mask - 1, 2 * mask]) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


basis = find(dirichlet == 0);
num_basis = size(basis, 1)
u_n = zeros(num_basis, 1);
k_val = ones(num_basis, 1); 

K = sparse(basis, [1:num_basis]', k_val, 2 * ver_num, num_basis);

% K = sparse(2 * ver_num, num_basis);
% 
% for i = 1 : num_basis
%    
%     K(basis(i), i) = 1;
%     
% end




% for iter = 1 : 100
%    
%     for i = 1 : ver_num
%        
%         if boundary_ver(i) == 0
%            
%             one_ring = find(exam_mesh.e(i, :) > 0);
%             
%             num = size(one_ring, 2);
%             
%             ox = 0;
%             oy = 0;
%             
%             for j = 1 : num
%                
%                 ox = ox + q(2 * one_ring(j) - 1);
%                 oy = oy + q(2 * one_ring(j));
%                 
%             end
%             
%             ox = ox / num;
%             oy = oy / num;
%             
%             q(2 * i - 1) = ox;
%             q(2 * i) = oy;
%             
%         end
%         
%     end
%     
%     plot_result(u_n, 0);
%     
% end









































global_cache_claim



tol_x_cnt = 0;
tol_f_cnt = 0;

tol_x = 1e-16;
tol_f = 1e-16;

stop_cnt = 5;

amips_s = 1;

cmax = 1;


if strcmp(en_type, 'arap')
   
    cmin = 0;
    
end

if strcmp(en_type, 'mips')
   
    cmin = 2;
    
end

if strcmp(en_type, 'iso')
   
    cmin = 4;
    
end

if strcmp(en_type, 'amips')
   
    cmin = exp(amips_s * 2);
    
end


if strcmp(en_type, 'conf')
   
    cmin = 1;
    
end


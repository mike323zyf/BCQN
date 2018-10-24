
global ver_num exam_mesh q_rest dirichlet en_type cmin cmax K q

colormap(jet);

en_type = 'iso';

load_input_mesh('../data/ani1');

q_n = reshape(exam_mesh.v(:, 1:2)', ver_num * 2, 1);

q_n(1:2:ver_num*2) = q_n(1:2:ver_num*2) * 0.2;

q_rest = q_n;

dirichlet = zeros(2 * ver_num, 1);

anchor_fix = find(bsxfun(@and, exam_mesh.v(:, 1) < -90, exam_mesh.v(:, 2) < 412));
anchor_move = find(exam_mesh.v(:, 1) > 1033);

ctr = mean(exam_mesh.v(anchor_move, 1:2));

theta = -90 * pi / 180;
rot = [cos(theta) sin(theta); -sin(theta) cos(theta)];
q_rest([2 * anchor_move - 1, 2 * anchor_move]) = bsxfun(@plus, bsxfun(@plus, bsxfun(@minus, exam_mesh.v(anchor_move, 1:2), ctr) * rot, ctr), [-600, -400]);
q_rest([2 * anchor_fix - 1, 2 * anchor_fix]) = exam_mesh.v(anchor_fix, 1:2);

mask = [anchor_fix; anchor_move];

dirichlet([2 * mask - 1, 2 * mask]) = 1;

q = q_rest;

basis = find(dirichlet == 0);
num_basis = size(basis, 1);
u_n = zeros(num_basis, 1);
k_val = ones(num_basis, 1); 

K = sparse(basis, [1:num_basis]', k_val, 2 * ver_num, num_basis);

global_cache_claim


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



global ver_num tet_mesh dirichlet en_type cmin cmax K q KK bb

colormap(jet);

en_type = 'iso';

load_input_mesh('../data/dancer1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q = reshape(tet_mesh.v(:, :)', ver_num * 3, 1);

% q = rand(ver_num * 3, 1);

anchor1 = find(tet_mesh.v(:, 2) < -144);
anchor2 = intersect(find(tet_mesh.v(:, 2) > 15), find(tet_mesh.v(:, 1) > 30));
anchor3 = intersect(find(tet_mesh.v(:, 2) > 40), find(tet_mesh.v(:, 1) < 10));

load('../data/dancer_init.mat');

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

soft_dirichlet = zeros(3 * ver_num, 1);

soft_mask = [anchor2; anchor3];

soft_dirichlet([3 * soft_mask - 2, 3 * soft_mask - 1, 3 * soft_mask]) = 1;


soft_basis = find(soft_dirichlet == 1);
soft_num_basis = size(soft_basis, 1);
soft_k_val = ones(soft_num_basis, 1); 

KK = sparse([1:soft_num_basis]', soft_basis, soft_k_val, soft_num_basis, 3 * ver_num);

bb = zyf(soft_basis);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

mask = [anchor1];

% q(3 * anchor2 - 2) = q(3 * anchor2 - 2) + 50;
% q(3 * anchor3 - 2) = q(3 * anchor3 - 2) - 40;
% q(3 * anchor3 - 1) = q(3 * anchor3 - 1) + 20;

dirichlet = zeros(3 * ver_num, 1);
dirichlet([3 * mask - 2, 3 * mask - 1, 3 * mask]) = 1;

basis = find(dirichlet == 0);
num_basis = size(basis, 1);
u_n = zeros(num_basis, 1);
k_val = ones(num_basis, 1); 

K = sparse(basis, [1:num_basis]', k_val, 3 * ver_num, num_basis);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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


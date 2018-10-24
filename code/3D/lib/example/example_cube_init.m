
global ver_num tet_mesh dirichlet en_type cmin cmax K q KK bb

colormap(jet);

en_type = 'gmr';

load_input_mesh('../data/soft_cube');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q = reshape(tet_mesh.v(:, :)', ver_num * 3, 1);

fix_anchor = find(tet_mesh.v(:, 2) < -0.19);
move_anchor = find(tet_mesh.v(:, 2) > 0.19);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

soft_dirichlet = zeros(3 * ver_num, 1);

soft_mask = [move_anchor];

soft_dirichlet([3 * soft_mask - 2, 3 * soft_mask - 1, 3 * soft_mask]) = 1;


soft_basis = find(soft_dirichlet == 1);
soft_num_basis = size(soft_basis, 1);
soft_k_val = ones(soft_num_basis, 1); 

KK = sparse([1:soft_num_basis]', soft_basis, soft_k_val, soft_num_basis, 3 * ver_num);

sq = q;
sq(3 * soft_mask - 1) = sq(3 * soft_mask - 1) + 1;

bb = sq(soft_basis);

% bb = q(3 * soft_mask - 1) + 1;
% bb = zyf(soft_basis);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

mask = [fix_anchor];

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



















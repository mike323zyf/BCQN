
global ver_num tet_mesh dirichlet en_type cmin cmax K q KK bb

colormap(jet);

en_type = 'iso';

load_input_mesh('../data/wrench');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q = reshape(tet_mesh.v(:, :)', ver_num * 3, 1);

fix_anchor = find(tet_mesh.v(:, 1) > 63);
move_anchor = find(tet_mesh.v(:, 1) < -70);

q_move = q;

theta = -100;
V = [q_move(3 * move_anchor - 2), q_move(3 * move_anchor - 1), q_move(3 * move_anchor)];
ctr = mean(V);
theta = theta * pi / 180;
rot = [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0; 0 0 1];

V = bsxfun(@plus,bsxfun(@minus,V,ctr)*rot,ctr);
q_move(3 * move_anchor - 2) = V(:, 1) + 100;
q_move(3 * move_anchor - 1) = V(:, 2) - 80;
q_move(3 * move_anchor - 0) = V(:, 3) - 140;

% q = q_move;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

soft_dirichlet = zeros(3 * ver_num, 1);

soft_mask = [move_anchor];

soft_dirichlet([3 * soft_mask - 2, 3 * soft_mask - 1, 3 * soft_mask]) = 1;


soft_basis = find(soft_dirichlet == 1);
soft_num_basis = size(soft_basis, 1);
soft_k_val = ones(soft_num_basis, 1); 

KK = sparse([1:soft_num_basis]', soft_basis, soft_k_val, soft_num_basis, 3 * ver_num);

sq = q_move;
% sq(3 * soft_mask - 1) = sq(3 * soft_mask - 1) + 1;

bb = sq(soft_basis);

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



















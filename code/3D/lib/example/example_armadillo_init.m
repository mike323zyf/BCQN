
global ver_num tet_mesh dirichlet en_type cmin cmax K q

colormap(jet);

en_type = 'iso';

load_input_mesh('../data/armadillo');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q = reshape(tet_mesh.v(:, :)', ver_num * 3, 1);

% q = rand(ver_num * 3, 1);

mask = [1, 174, 3962, 8283, 10293, 11570];


q(174 * 3 - 2) = q(174 * 3 - 2) + 32;
q(174 * 3) = q(174 * 3) - 16; 

q(3962 * 3 - 2) = q(3962 * 3 - 2) + 8;
q(3962 * 3 - 1) = q(3962 * 3 - 1) - 3;
q(3962 * 3) = q(3962 * 3) + 18; 

q(8283 * 3 - 2) = q(8283 * 3 - 2) + 40;
q(8283 * 3 - 0) = q(8283 * 3 - 0) + 40;

q(10293 * 3 - 2) = q(10293 * 3 - 2) + 20;
q(10293 * 3 - 0) = q(10293 * 3 - 0) - 20;

q(11570 * 3 - 2) = q(11570 * 3 - 2) + 12;
q(11570 * 3 - 1) = q(11570 * 3 - 1) + 30;

dirichlet = zeros(3 * ver_num, 1);
dirichlet([3 * mask - 2, 3 * mask - 1, 3 * mask]) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


basis = find(dirichlet == 0);
num_basis = size(basis, 1);
u_n = zeros(num_basis, 1);
k_val = ones(num_basis, 1); 

K = sparse(basis, [1:num_basis]', k_val, 3 * ver_num, num_basis);


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


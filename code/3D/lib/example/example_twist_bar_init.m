
global ver_num tet_mesh dirichlet en_type cmin cmax K q

colormap(jet);

en_type = 'gmr';

load_input_mesh('../data/rest');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q = reshape(tet_mesh.v(:, :)', ver_num * 3, 1);

anchor1 = find(q > 0.24) / 3;
anchor2 = find(q < -0.24) / 3;

mask = [anchor1; anchor2];

theta = 42;
V = [q(3 * anchor1 - 2), q(3 * anchor1 - 1), q(3 * anchor1)];
ctr = mean(V);
theta = theta * pi / 180;
rot = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];

V = bsxfun(@plus,bsxfun(@minus,V,ctr)*rot,ctr);
q(3 * anchor1 - 2) = V(:, 1);
q(3 * anchor1 - 1) = V(:, 2);
q(3 * anchor1 - 0) = V(:, 3);


theta = -44;
V = [q(3 * anchor2 - 2), q(3 * anchor2 - 1), q(3 * anchor2)];
ctr = mean(V);
theta = theta * pi / 180;
rot = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];

V = bsxfun(@plus,bsxfun(@minus,V,ctr)*rot,ctr);
q(3 * anchor2 - 2) = V(:, 1);
q(3 * anchor2 - 1) = V(:, 2);
q(3 * anchor2 - 0) = V(:, 3);

q(anchor1 * 3) = q(anchor1 * 3) + 0.0;
q(anchor2 * 3) = q(anchor2 * 3) - 0.1;

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


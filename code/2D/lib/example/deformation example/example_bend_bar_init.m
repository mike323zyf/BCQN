global ver_num tri_num obj_tri exam_mesh dirichlet en_type cmin cmax K q

colormap(jet);

% en_type = 'gmr'
% en_type = 'mips'
en_type = 'iso'
% en_type = 'arap'
% en_type = 'conf' % don't do conf for defo

exam_mesh = struct('v', [], 'n', [], 'u', [], 'f', [], 'e', [], 'bounds', []);

n = 100;

n_x = n;
n_y = ceil(n/20.0);

[x, y] = meshgrid(1:n_x, 1:n_y);
x = x(:); y = y(:);
V = [x, y];
F = delaunay(V);

exam_mesh.v = V;
exam_mesh.f = F;

ver_num = size(exam_mesh.v, 1);

tri_num = size(exam_mesh.f, 1);

obj_tri = exam_mesh.f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dirichlet = zeros(2 * ver_num, 1);

anchors_fix = find(x<2);
anchors_move = find(x==n_x);

mask = [anchors_fix; anchors_move];

dirichlet([2 * mask - 1, 2 * mask]) = 1;

for i = 1 : n
   anchors = find(x == i);
   theta = (i - 1) * pi / (n - 1);
   rot = [cos(theta) sin(theta); -sin(theta) cos(theta)];
   q_n([2 * anchors - 1, 2 * anchors]) = bsxfun(@minus, V(anchors_fix, :), [1, 1 + n / 20]) * rot;
end

q = q_n';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basis = find(dirichlet == 0);
num_basis = size(basis, 1);
u_n = zeros(num_basis, 1);
k_val = ones(num_basis, 1); 

K = sparse(basis, [1:num_basis]', k_val, 2 * ver_num, num_basis);


global_cache_claim


cmax = 0.5;


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

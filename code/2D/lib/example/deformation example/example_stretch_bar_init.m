global ver_num tri_num obj_tri exam_mesh dirichlet en_type cmin cmax K q c1 c2

colormap(jet);

en_type = 'gmr'
% en_type = 'mips'
% en_type = 'iso'
%en_type = 'arap'
% en_type = 'conf' % don't do conf for defo

exam_mesh = struct('v', [], 'n', [], 'u', [], 'f', [], 'e', [], 'bounds', []);

%n=100;
n=50;
n_y = n;
n_x = ceil(n/2.0);

[x, y] = meshgrid(1:n_x, 1:n_y);
x = x(:); y = y(:);

cx = mean(x);
cy = mean(y);
x = (x - cx) / (n_x - 1);
y = (y - cy) / (n_y - 1) * 2;
miny = min(y);
maxy = max(y);

V = [x, y];
F = delaunay(V);

exam_mesh.v = V;
exam_mesh.f = F;

ver_num = size(exam_mesh.v, 1);

tri_num = size(exam_mesh.f, 1);

obj_tri = exam_mesh.f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q = reshape(exam_mesh.v(:, 1:2)', ver_num * 2, 1);

q_rest = q;

for i = 1 : ver_num
   
    q(2 * i) = q(2 * i) * 0.5 - 0.5;
    
end

dirichlet = zeros(2 * ver_num, 1);

anchors_fix = find(y==miny);
anchors_move = find(y==maxy);

mask = [anchors_fix; anchors_move];

dirichlet([2 * mask - 1, 2 * mask]) = 1;

q_rest(2 * anchors_move) = q_rest(2 * anchors_move) + 3;

q(find(dirichlet)) = q_rest(find(dirichlet));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basis = find(dirichlet == 0);
num_basis = size(basis, 1);
u_n = zeros(num_basis, 1);
k_val = ones(num_basis, 1); 

K = sparse(basis, [1:num_basis]', k_val, 2 * ver_num, num_basis);


global_cache_claim


cmax = 2;


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

if strcmp(en_type, 'gmr')
   
    cmin = -c1 - 2 * c2;
    
end

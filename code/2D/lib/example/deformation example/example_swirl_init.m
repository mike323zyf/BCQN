global ver_num tri_num obj_tri exam_mesh q_rest dirichlet en_type cmin cmax K q


% en_type = 'gmr'
% en_type = 'mips'
en_type = 'arap'
% en_type = 'iso'
%en_type = 'conf' % don't do conf for defo

exam_mesh = struct('v', [], 'n', [], 'u', [], 'f', [], 'e', [], 'bounds', []);


total_theta = 8 * pi;
delta_radius = 1.5;
total_radius = 1.0;

% n = 6400;
% n = 300;
n = 100;
% n = 3400;

n_x = n;
n_y = ceil(n/20.0);

step_theta = total_theta / (n_x - 1);
step_radius = total_radius / (n_y - 1);

[x, y] = meshgrid(1:n_x, 1:n_y);
x = x(:); y = y(:);
V = [x, y];
F = delaunay(V);

num = size(x, 1);

index = num * 2 * x + y;

tx = x;
ty = y;

fixed = find(x == 1);

mask = [fixed(1)];


for i = 1 : n_x
    
    theta = (i - 1) * step_theta + pi / 2;
    
    for j = 1 : n_y
        
        current_pt = find(index == (num * 2 * i + j));
        
        radius = (delta_radius * (1 + (i - 1) / (n_x - 1)) + (j - 1) * step_radius) * theta;
        
        tx(current_pt) = radius * cos(theta);
        ty(current_pt) = radius * sin(theta);
        
    end
    
end

x = tx;
y = ty;

V = [x, y];

exam_mesh.v = V;
exam_mesh.f = F;

ver_num = size(exam_mesh.v, 1);

tri_num = size(exam_mesh.f, 1);

obj_tri = exam_mesh.f;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

dirichlet = zeros(2 * ver_num, 1);

dirichlet([2 * mask - 1, 2 * mask]) = 1;


for i = 1 : n_x
    
    theta = (1 - i) * step_theta - pi / 2;
    
    for j = 1 : n_y
        
        current_pt = find(index == (num * 2 * i + (n_y + 1 - j)));
        
        radius = (delta_radius * (1 + (i - 1) / (n_x - 1)) + (j - 1) * step_radius) * abs(theta);
        
        tx(current_pt) = radius * cos(theta);
        ty(current_pt) = radius * sin(theta);
        
    end
    
end

V = [tx, ty];

q = reshape(V', ver_num * 2, 1);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

basis = find(dirichlet == 0);
num_basis = size(basis, 1);
u_n = zeros(num_basis, 1);
k_val = ones(num_basis, 1); 

K = sparse(basis, [1:num_basis]', k_val, 2 * ver_num, num_basis);

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

global_cache_claim

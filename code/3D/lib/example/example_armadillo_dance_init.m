

global ver_num tet_mesh dirichlet en_type cmin cmax K q KK bb

colormap(jet);

en_type = 'gmr';

load_input_mesh('../data/armadillo');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q = reshape(tet_mesh.v(:, :)', ver_num * 3, 1);

fix_anchor = intersect(find(tet_mesh.v(:, 1) > 0), find(tet_mesh.v(:, 2) < -40));
move_anchor1 = intersect(find(tet_mesh.v(:, 1) < 0), find(tet_mesh.v(:, 2) < -40));
move_anchor2 = intersect(intersect(find(tet_mesh.v(:, 1) < -40), find(tet_mesh.v(:, 2) > 30)), find(tet_mesh.v(:, 3) < -30));
move_anchor3 = intersect(intersect(find(tet_mesh.v(:, 1) > 40), find(tet_mesh.v(:, 2) > 30)), find(tet_mesh.v(:, 3) < -20));
% move_anchor4 = intersect(intersect(intersect(find(tet_mesh.v(:, 1) > -30), find(tet_mesh.v(:, 1) < 30)), find(tet_mesh.v(:, 2) > 30)), find(tet_mesh.v(:, 3) < -18));

q_move = q;

theta = -80;
V = [q_move(3 * move_anchor1 - 2), q_move(3 * move_anchor1 - 1), q_move(3 * move_anchor1)];
ctr = mean(V);
theta = theta * pi / 180;
rot = [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0; 0 0 1];
% rot = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];

V = bsxfun(@plus,bsxfun(@minus,V,ctr)*rot,ctr);
q_move(3 * move_anchor1 - 2) = V(:, 1) + 30;
q_move(3 * move_anchor1 - 1) = V(:, 2) + 30;
q_move(3 * move_anchor1 - 0) = V(:, 3);



theta = -80;
V = [q_move(3 * move_anchor2 - 2), q_move(3 * move_anchor2 - 1), q_move(3 * move_anchor2)];
ctr = mean(V);
theta = theta * pi / 180;
% rot = [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0; 0 0 1];
% rot = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
rot = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];

V = bsxfun(@plus,bsxfun(@minus,V,ctr)*rot,ctr);
q_move(3 * move_anchor2 - 2) = V(:, 1) - 40;
q_move(3 * move_anchor2 - 1) = V(:, 2);
q_move(3 * move_anchor2 - 0) = V(:, 3) + 50;



theta = -80;
V = [q_move(3 * move_anchor3 - 2), q_move(3 * move_anchor3 - 1), q_move(3 * move_anchor3)];
ctr = mean(V);
theta = theta * pi / 180;
% rot = [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0; 0 0 1];
rot = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];

V = bsxfun(@plus,bsxfun(@minus,V,ctr)*rot,ctr);
q_move(3 * move_anchor3 - 2) = V(:, 1) - 30;
q_move(3 * move_anchor3 - 1) = V(:, 2) + 60;
q_move(3 * move_anchor3 - 0) = V(:, 3) + 40;



% theta = -80;
% V = [q_move(3 * move_anchor4 - 2), q_move(3 * move_anchor4 - 1), q_move(3 * move_anchor4)];
% ctr = mean(V);
% theta = theta * pi / 180;
% % rot = [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0; 0 0 1];
% rot = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
% 
% V = bsxfun(@plus,bsxfun(@minus,V,ctr)*rot,ctr);
% q_move(3 * move_anchor4 - 2) = V(:, 1);
% q_move(3 * move_anchor4 - 1) = V(:, 2) + 30;
% q_move(3 * move_anchor4 - 0) = V(:, 3) + 35;


% q = q_move;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

soft_dirichlet = zeros(3 * ver_num, 1);

soft_mask = [move_anchor1; move_anchor2; move_anchor3];

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


















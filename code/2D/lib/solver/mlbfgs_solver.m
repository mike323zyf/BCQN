function [ ] = mlbfgs_solver( un, is_uv_mesh )
%PPAP_SOLVER Summary of this function goes here
%   Detailed explanation goes here

global order total_volume

u = un;
p = 0 * u;

num = size(u, 1);

m = 5;

ss = zeros(num, m);
yy = zeros(num, m);

A = get_precomputed_info(is_uv_mesh);

scale = condest(A);
scale = scale^(1.0/3);

multiplier = 1;

order = symamd(A);
As = A(order, order);
Ac = chol(As);

plot_result( u, 0 )

grad = grad_function(u);

time = 0;
iteration = -1;

for i = 0 : 60000
    iteration = i;
    
%     save_mesh(u, 2, 3, i);

%     energy_value(u)

    plot_result( u, i )
    
    tic;
    grads = grad; 

    p(order) = lbfgs_update(-1.0 * grads, ss(:, max(m + 1 - i, 1) : m), yy(:, max(m + 1 - i, 1) : m), Ac);
    
    
    %%%
    
    p_mod = lcp_solve(p, 20, 0.5);
    
    p_new = p + p_mod;
    
    check = grads' * (p_new);
    
    if check <= 0
        
        p = p_new;
        
    end

    %%%
    
    
    un = line_check_search(p, u, grads); 

    if stop_check(un, u, grads)
        
        break;
        
    end   
    
    grad = grad_function(un);
    sk = un - u;
    yk = grad - grads;
    
    ss(:, 1 : m - 1) = ss(:, 2 : m);
    ss(:, m) = sk;
    yy(:, 1 : m - 1) = yy(:, 2 : m);
    
    zk = A * sk;
    
    %ble = (yk' * zk) / total_volume * scale * multiplier;
    ble = (yk' * zk) / total_volume * scale;
    %ble = (yk' * zk) / total_volume  * 1/norm(yk); 
    
%     disp('y*L*p :'); disp(yk' * zk);
%     disp('y*p :'); disp(sk'*yk);
%     disp('ble :'); disp(ble);
    
    ble = max(min(1.0, ble), 0); % Blend
%     ble = 1 % Sobolev
%     ble = 0 % BFGS
    
    
    yy(:, m) = ble * zk + (1 - ble) * yk;
    
    u = un;
    
    sub_time = toc;
    time = time + sub_time;
   
end

% stop_energy = energy_value(un)
% 
% stop_gradient = norm(grads) / total_volume
% 
% iteration
% 
% time

end


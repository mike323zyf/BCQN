function [ ] = mlbfgs_lcp_solver( un )
%PPAP_SOLVER Summary of this function goes here
%   Detailed explanation goes here

global ver_num order total_volume soft_param

u = un;
p = 0 * u;

num = size(u, 1);

m = 5;

ss = zeros(num, m);
yy = zeros(num, m);

A = get_precomputed_info();

scale = condest(A);
scale = scale^(1.0/3)
compare = normest(A)
% scale = 1;

multiplier = 1;

order = symamd(A);
As = A(order, order);
Ac = chol(As);
AcT = Ac';

% plot_result( u, 0 )

grad = grad_function(u);

time = 0;
iteration = -1;

for i = 0 : 1200
    iteration = i
    
%     save_mesh(u, i);
    
    tic;
    grads = grad; 

    p = opt_dual_update(grads, ss(:, max(m + 1 - i, 1) : m), yy(:, max(m + 1 - i, 1) : m), Ac, AcT);
    
    un = line_check_search(p, u, grads); 
    
%     grad_norm = norm(grads)
    
%     if grad_norm < 1e-2
%         
%         break;
%         
%     end

    if stop_check(un, u)
        
        break;
        
    end   
    
    grad = grad_function(un);
    sk = un - u;
    yk = grad - grads;
    
    ss(:, 1 : m - 1) = ss(:, 2 : m);
    ss(:, m) = sk;
    yy(:, 1 : m - 1) = yy(:, 2 : m);
    
    zk = A * sk;
    
    ble = (yk' * zk) / total_volume^(4.0 / 3.0) * scale * multiplier
    ble = max(min(1.0, ble), 0);
    
    yy(:, m) = ble * zk + (1 - ble) * yk;
    
    u = un;
    
    sub_time = toc;
    time = time + sub_time;
    
    plot_result( u, i )
   
end

plot_result( u, i )

stop_energy = energy_value(un)

stop_gradient = norm(grads) / total_volume

iteration

time

end


function [ ] = mlbfgs_lcp_split_solver( un )
%PPAP_SOLVER Summary of this function goes here
%   Detailed explanation goes here

global order total_volume free_num

u = un;
p = 0 * u;
zk = p;

num = size(u, 1);

m = 5;

ss = zeros(num, m);
yy = zeros(num, m);

A = get_precomputed_split_info();

scale = normest(A);

Ac = chol(A);
AcT = Ac';

grad = grad_function(u);

for i = 0 : 100000

    i

    plot_result( u, i )
    
    grads = grad; 
  
    p = opt_dual_update(grads, ss(:, max(m + 1 - i, 1) : m), yy(:, max(m + 1 - i, 1) : m), Ac, AcT);

    un = line_check_search(p, u, grads); 
 
    if stop_check(un, u, grads)
        
        break;
        
    end   

    grad = grad_function(un);

    sk = un - u;
    yk = grad - grads;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ss(:, 1 : m - 1) = ss(:, 2 : m);
    ss(:, m) = sk(order);
    yy(:, 1 : m - 1) = yy(:, 2 : m);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    zk(order(1:free_num)) = A * sk(order([1:free_num]));
    zk(order(free_num+1:2*free_num)) = A * sk(order(free_num+1:2*free_num));
    zk(order(2*free_num+1:3*free_num)) = A * sk(order(2*free_num+1:3*free_num));
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ble = (yk' * zk) / total_volume^(4.0 / 3.0) * scale;
    ble = max(min(1.0, ble), 0);
   
    tt = ble * zk + (1 - ble) * yk;
    yy(:, m) = tt(order);
    
    u = un;
    
end

end


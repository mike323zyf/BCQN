function [ ] = qp_solver( un, is_uv_mesh )
%BFGS_SOLVER Summary of this function goes here
%   Detailed explanation goes here

global ver_num order total_volume

tic;

p = un * 0;
zyf = p;

A = get_precomputed_info(is_uv_mesh);

time1 = toc


tic;

order = symamd(A);
As = A(order, order);
Ac = chol(As);
AcT = Ac';

time2 = toc

tic;

[row, col, val] = find(As);

time3 = toc

nnz_amg = size(val, 1)
dof_amg = size(As, 1)

maxNumCompThreads(8);

tic;

lapamg_solve_mex(nnz_amg, dof_amg, col, row, 1.0 * val, 0, rand(dof_amg, 1));

time4 = toc

for i = 1 : 2
    i

    grads = grad_function(un); 

    tic;
    p(order) = -1.0 * lapamg_solve_mex(nnz_amg, dof_amg, col, row, 1.0 * val, 1, grads(order));
    toc
    
    tic;
    zyf(order) = -1.0 * (Ac \ (AcT \ grads(order))); 
    toc
    
    norm(zyf - p) / norm(zyf)
    
    us = line_check_search(p, un, grads); 

    if stop_check(us, un, grads)
        
        break;
        
    end
       
    un = us;
    
    
    
    plot_result( un, i )
    
end

end


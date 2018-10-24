function [ ] = aqp_solver( un, is_uv_mesh )
%BFGS_SOLVER Summary of this function goes here
%   Detailed explanation goes here

global ver_num order total_volume

u_n0 = un;
u_n1 = un;
u_n2 = un;
p = un * 0;

% theta = 0.8182;
% theta = 0.9387;
theta = 0.6;

num = size(u_n1, 1);

A = get_precomputed_info(is_uv_mesh);

order = symamd(A);
As = A(order, order);
Ac = chol(As);
AcT = Ac';

% plot_result( u_n0, 0 )

time = 0;

for i = 0 : 1000
    i
    
%     if mod(i, 4) == 0
        plot_result( u_n1, i )
%     end
    tic;
    y_n = get_feasible_acceleration(u_n1 - u_n0, u_n1, theta); 
    
    grads = grad_function(y_n); 
    
    p(order) = -1.0 * (Ac \ (AcT \ grads(order))); 
    
    u_n2 = line_check_search(p, y_n, grads); 

    if stop_check(u_n2, u_n1, grads)
        
        break;
        
    end
       
    u_n0 = u_n1;
    u_n1 = u_n2;
    
    time = time + toc;
    
    if mod(i, 100) == 0
%         
        
%         
    end
    
     

    
end

% plot_result( u_n1, i )

i

time

% plot_result( u_n1, i )

end


function [ ] = aqp_solver( un )
%BFGS_SOLVER Summary of this function goes here
%   Detailed explanation goes here

global ver_num order

u_n0 = un;
u_n1 = un;
u_n2 = un;
p = un * 0;

theta = 0.8182;
% theta = 0.9387;
% theta = 0.0;

A = get_precomputed_info();

order = symamd(A);
As = A(order, order);
Ac = chol(As);
AcT = Ac';

% plot_result( u_n0, 0 )
% save_mesh(u_n1, 0);

time = 0;

for i = 0 : 1200
    
    i
    
%     save_mesh(u_n1, i);
    tic;
    y_n = get_feasible_acceleration(u_n1 - u_n0, u_n1, theta);
    
    grads = grad_function(y_n); 
        
    p(order) = -1.0 * (Ac \ (AcT \ grads(order)));  
    
    u_n2 = line_check_search(p, y_n, grads);
    
%     grad_norm = norm(grads)

%     check = energy_value(u_n2)
    
%     if check < 3.9940e+05
%         
%         break;
%         
%     end
    
    if stop_check(u_n2, u_n1, grads)
        
        break;
        
    end
       
    u_n0 = u_n1;
    u_n1 = u_n2;
    
    sub_time = toc;
    time = time + sub_time
    
%     if mod(i, 10) == 0
%         plot_result( u_n1, i )
%     end
   
end

time

% plot_result( u_n1, i )

end


function [ ] = aqp_split_solver( un, is_uv_mesh )
%BFGS_SOLVER Summary of this function goes here
%   Detailed explanation goes here

global order order

u_n0 = un;
u_n1 = un;
u_n2 = un;
p = un * 0;

% theta = 0.8182;
theta = 0.9387;
% theta = 0.5;

A = get_precomputed_split_info(is_uv_mesh);

Ac = chol(A);
AcT = Ac';

for i = 0 : 250
    
    i
    
    plot_result( u_n1, i)
   
    y_n = get_feasible_acceleration(u_n1 - u_n0, u_n1, theta); 
    
    grads = grad_function(y_n); 
    
    p(order) = -1.0 * split_solve(Ac, AcT, grads(order));
    
    u_n2 = line_check_search(p, y_n, grads); 

    if stop_check(u_n2, u_n1, grads)
        
        break;
        
    end
       
    u_n0 = u_n1;
    u_n1 = u_n2;
    
end

end


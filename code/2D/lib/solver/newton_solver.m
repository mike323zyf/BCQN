function [ ] = newton_solver( un, is_uv_mesh )
%NEWTON_SOLVER Summary of this function goes here
%   Detailed explanation goes here

u = un;
p = 0 * u;

L = get_precomputed_info(is_uv_mesh);

order = get_ordering(u);

'newton'

for i = 0 : 180
    
    i
    
    plot_result( u, i )

    [grad, H] = grad_hessian_function(u, 0);

    p(order) = -1.0 * H(order, order) \ grad(order);

    un = line_check_search(p, u, grad);
    
    if stop_check(un, u, grad)
        
        break;
        
    end  
    
    u = un;
    
end

end


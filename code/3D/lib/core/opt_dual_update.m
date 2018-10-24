function [ p ] = opt_dual_update( g, s, y, Ac, AcT )
%DUAL_PPAP Summary of this function goes here
%   Detailed explanation goes here

global order

p = g * 0;

p(order) = lbfgs_split_update(-1.0 * g, s, y, Ac, AcT);

p_mod = lcp_solve(p, 20, 0.5);

p_new = p + p_mod;

check = g' * (p_new);

if check <= 0
    
    p = p_new;
    
end

end


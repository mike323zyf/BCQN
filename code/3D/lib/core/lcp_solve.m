function [ p_mod ] = lcp_solve( p, iter, lcp_a )
%LCP_SOLVE Summary of this function goes here
%   Detailed explanation goes here

global J_u J_ui JT_u JT_ui JV_u JTV_u wu bu USE_TBB

JTJ_info = [J_ui(1), J_ui(2), JT_ui(1), JT_ui(2)];

if USE_TBB == 1
    
    maxNumCompThreads(8);
    
    p_mod = lcp_solve_tbb_mex(J_u, JT_u, JV_u, JTV_u, JTJ_info, p, wu, bu, iter, lcp_a);
    
else
    
    p_mod = lcp_solve_eigen_mex(J_u, JT_u, JV_u, JTV_u, JTJ_info, p, wu, bu, iter, lcp_a);
    
end

end


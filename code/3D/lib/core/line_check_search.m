function [ output ] = line_check_search( p, u, grad )
%LINE_CHECK_SEARCH Summary of this function goes here
%   Detailed explanation goes here

if get_energy_type() > 0
    
    fk = energy_value(u);
    dfk = grad;

    ls_alpha = 0.2;
    ls_beta = 0.5;

    balp = get_feasible_step_size(p, u);
    alp = min(1.0, balp);

    un = u + alp * p;
    lhs = energy_value(un);
    rhs = fk + ls_alpha * alp * dot(dfk, p);

    while lhs > rhs
    
        alp = alp * ls_beta;
    
        un = u + alp * p;
        lhs = energy_value(un);
        rhs = fk + ls_alpha * alp * dot(dfk, p);
        
    end
    
    output = un;
    
else
    
    output = u + p;
    
end

end


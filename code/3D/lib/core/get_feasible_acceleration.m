function [ y_n ] = get_feasible_acceleration( p, u, theta )
%GET_FEASIBLE_ACCELERATION Summary of this function goes here
%   Detailed explanation goes here

global K q tet_mesh ver_num

if get_energy_type() > 0
    
    q_target = K * u + q;

    p_direction = K * p;

    qq = reshape(q_target, 3, ver_num)';
    pp = reshape(p_direction, 3, ver_num)';

    q_target = reshape(qq, 3 * ver_num, 1);
    p_direction = reshape(pp, 3 * ver_num, 1);

    tolerance = 1e-8;

    steps = get_feasible_steps_mex(tet_mesh.t, q_target, p_direction, tolerance);

    output = min(steps);

    output = 0.5 * output;

    theta = min(output, theta);
    
end

y_n = u + theta * p;

end


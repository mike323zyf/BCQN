function [ output ] = get_feasible_step_size( p, u )
%GET_FEASIBLE_STEP_SIZE Summary of this function goes here
%   Detailed explanation goes here

global K q exam_mesh ver_num

q_target = K * u + q;

p_direction = K * p;

qq = reshape(q_target, 2, ver_num)';
pp = reshape(p_direction, 2, ver_num)';

q_target = reshape(qq, 2 * ver_num, 1);
p_direction = reshape(pp, 2 * ver_num, 1);

tolerance = 1e-5;

steps = get_feasible_steps_mex(exam_mesh.f, q_target, p_direction, tolerance);

output = min(steps);

output = 0.5 * output;

end


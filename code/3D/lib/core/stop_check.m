function [ output ] = stop_check( un2, un1, grad )
%STOP_CHECK Summary of this function goes here
%   Detailed explanation goes here

global tol_x_cnt tol_f_cnt stop_cnt tol_x tol_f perimeter_norm

output = 0;

if norm(un2 - un1) < tol_x * (1 + norm(un1))
    tol_x_cnt = tol_x_cnt + 1; 
else
    tol_x_cnt = 0;
end

if tol_x_cnt >= stop_cnt
    output = 1; 
end

fn2 = energy_value(un2);
fn1 = energy_value(un1);

if abs(fn1 - fn2) < tol_f * (1 + abs(fn1))
    tol_f_cnt = tol_f_cnt + 1;
else
    tol_f_cnt = 0;
end

if tol_f_cnt >= stop_cnt
    output = 1; 
end

if norm(grad) < 1e-3 * perimeter_norm
    output = 1;
    
    stopped_type = 3    
end


end


function [ lhs ] = split_solve( A, AT, rhs )
%SPLIT_SOLVE Summary of this function goes here
%   Detailed explanation goes here

global free_num

% 'split solve'

rhs1 = rhs([1:free_num]);
rhs2 = rhs([free_num+1:2*free_num]);
rhs3 = rhs([2*free_num+1:3*free_num]);

lhs_c = [rhs1, rhs2, rhs3];

lhs_c(:, :) = A \ (AT \ [rhs1, rhs2, rhs3]);

lhs = rhs;

lhs([1:free_num]) = lhs_c(:, 1);
lhs([free_num+1:2*free_num]) = lhs_c(:, 2);
lhs([2*free_num+1:3*free_num]) = lhs_c(:, 3);

end


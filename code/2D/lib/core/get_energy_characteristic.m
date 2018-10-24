function [ output ] = get_energy_characteristic( )
%GET_ENERGY_CHARACTERISTIC Summary of this function goes here
%   Detailed explanation goes here

global en_type c1 c2 d1

hessian = zeros(2, 2);

if strcmp(en_type, 'arap')
   
    hessian(1, 1) = 2;
    hessian(2, 2) = 2;
        
    output = norm(hessian);
    
end

if strcmp(en_type, 'mips')
   
    hessian(1, 1) = 2;
    hessian(2, 2) = 2;
    hessian(1, 2) = -2;
    hessian(2, 1) = -2;
        
    output = norm(hessian);
    
end

if strcmp(en_type, 'iso')
   
    hessian(1, 1) = 8;
    hessian(2, 2) = 8;
        
    output = norm(hessian);
    
end

if strcmp(en_type, 'amips')
   
    output = 3;
    
end

temp = 25;

if strcmp(en_type, 'conf')
   
    hessian(1, 1) = 2 * temp;
    hessian(2, 2) = 6 * temp;
    hessian(1, 2) = -4 * temp;
    hessian(2, 1) = -4 * temp;
        
    output = norm(hessian);
    
end


if strcmp(en_type, 'gmr')
   
    hessian(1, 1) = 2 * (c1 + d1);
    hessian(2, 2) = 2 * (c1 + d1);
    hessian(1, 2) = 2 * (d1 - c1);
    hessian(2, 1) = 2 * (d1 - c1);
        
    output = norm(hessian);
    
end

end


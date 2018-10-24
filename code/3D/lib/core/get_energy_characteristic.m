function [ output ] = get_energy_characteristic( )
%GET_ENERGY_CHARACTERISTIC Summary of this function goes here
%   Detailed explanation goes here

global en_type c1 c2 d1

hessian = zeros(3, 3);

if strcmp(en_type, 'arap')
   
    hessian(1, 1) = 2;
    hessian(2, 2) = 2;
    hessian(3, 3) = 2;
        
    output = norm(hessian);
    
end

if strcmp(en_type, 'mips')
   
    
    
end

if strcmp(en_type, 'iso')
   
    hessian(1, 1) = 8;
    hessian(2, 2) = 8;
    hessian(3, 3) = 8;
        
    output = norm(hessian);
    
end

if strcmp(en_type, 'amips')
   
    
    
end


if strcmp(en_type, 'conf')
   
    
    
end


if strcmp(en_type, 'gmr')
   
    hessian(1, 1) = 8 * (c1 + c2) / 3 + 2 * d1;
    hessian(1, 2) = -4 * (c1 + c2) / 3 + 2 * d1;
    hessian(2, 1) = -4 * (c1 + c2) / 3 + 2 * d1;
    hessian(1, 3) = -4 * (c1 + c2) / 3 + 2 * d1;
    hessian(3, 1) = -4 * (c1 + c2) / 3 + 2 * d1;
    hessian(2, 2) = 8 * (c1 + c2) / 3 + 2 * d1;
    hessian(2, 3) = -4 * (c1 + c2) / 3 + 2 * d1;
    hessian(3, 2) = -4 * (c1 + c2) / 3 + 2 * d1;
    hessian(3, 3) = 8 * (c1 + c2) / 3 + 2 * d1;
        
    output = norm(hessian);
    
end

end


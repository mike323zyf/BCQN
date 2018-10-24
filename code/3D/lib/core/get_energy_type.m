function [ output ] = get_energy_type()
%GET_ENERGY_TYPE Summary of this function goes here
%   Detailed explanation goes here

global en_type

if strcmp(en_type, 'arap')
   
    output = 0;
    
end

if strcmp(en_type, 'mips')
   
    output = 1;
    
end

if strcmp(en_type, 'iso')
   
    output = 2;
    
end

if strcmp(en_type, 'amips')
   
    output = 3;
    
end


if strcmp(en_type, 'conf')
   
    output = 4;
    
end


if strcmp(en_type, 'gmr')
   
    output = 5;
    
end

end


function [d] = lbfgs_update( gg, ss, yy, Ac, AcT )

global order

s = ss(order, :);
y = yy(order, :);
g = gg(order, :);

[p, k] = size(s);

for i = 1 : k
    
    iro = y(:, i)' * s(:, i);
    
    if abs(iro) < 1e-8
        
        iro = 1e-8;
        
    end
    
    ro(i, 1) = 1 / iro;
end

q = zeros(p, k + 1);
r = zeros(p, k + 1);
al =zeros(k, 1);
be =zeros(k, 1);

q(:, k + 1) = g;

for i = k : -1 : 1
    al(i) = ro(i) * s(:, i)' * q(:, i + 1);
    q(:, i) = q(:, i + 1) - al(i) * y(:, i);
end

r(:, 1) = Ac \ (AcT \ q(:, 1));

for i = 1 : k
    be(i) = ro(i) * y(:, i)' * r(:, i);
    r(:, i + 1) = r(:, i) + s(:, i) * (al(i) - be(i));
end

d = r(:, k + 1);

end

function [d] = lbfgs_split_update( gg, ss, yy, Ac, AcT )

global order

[p, k] = size(ss);

for i = 1 : k
    
    iro = yy(:, i)' * ss(:, i);
    
    if abs(iro) < 1e-8
        
        iro = 1e-8;
        
    end
    
    ro(i, 1) = 1 / iro;
end

q = zeros(p, k + 1);
r = zeros(p, k + 1);
al = zeros(k, 1);
be = zeros(k, 1);

q(:, k + 1) = gg(order);

for i = k : -1 : 1
    al(i) = ro(i) * (ss(:, i)' * q(:, i + 1));
    q(:, i) = q(:, i + 1) - al(i) * yy(:, i);
end

r(:, 1) = split_solve(Ac, AcT, q(:, 1));

for i = 1 : k
    be(i) = ro(i) * yy(:, i)' * r(:, i);
    r(:, i + 1) = r(:, i) + ss(:, i) * (al(i) - be(i));
end

d = r(:, k + 1);

end

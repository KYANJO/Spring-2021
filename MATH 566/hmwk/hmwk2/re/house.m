function [V,R] = house(A)

[m,n] = size(A);

V = zeros(m,n);

for k = 1:n
    
    x = A(k:m,k);
    l = length(x);
    vk = sign(x(1))*norm(x,2)*eye(l,1) + x;
    vk = vk/norm(vk,2);
    V(k:m,k) = vk;
    
    A(k:m,k:n) = A(k:m,k:n) - 2*vk*(vk'*A(k:m,k:n));
    
    
end

R = A;
end
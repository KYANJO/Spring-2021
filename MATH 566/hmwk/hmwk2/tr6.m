
%m equally spaced points over [0,1]
m = 5; n=4;

% Vandermonde matrix B
t = zeros(m,n);
for i = 1:n
    for j = 1:m
        t(j,i) = ((j-1)/(m-1))^(n-i);
    end
end

%fliping the vandermonde matrix t to form A
A = fliplr(t);

[V,R] = house(A);

Q = house2q(V)


function Q = house2q(V)

[m,n] = size(V);

x = eye(m);

for k = n:1
    x(k:m,:) = x(k:m,:) - 2*V(k:m,k)*(V(k:m,k)'*x(k:m,:));
end
Q = x;
end

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
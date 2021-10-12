clear all;
%m equally spaced points over [0,1]
m = 4; n=4;

% Vandermonde matrix B
t = zeros(m,n);
for i = 1:n
    for j = 1:m
        t(j,i) = ((j-1)/(m-1))^(n-i);
    end
end

%fliping the vandermonde matrix t to form A
A = fliplr(t);

Q = house2q(V);

function Q = house2q(V)

[m,n] = size(V);  

x = eye(m);

for k = n:1
    x(k:m,:) = x(k:m,:) - 2*V(k:m,k)*(V(k:m,k)'*x(k:m,:));
end
Q = x;
end
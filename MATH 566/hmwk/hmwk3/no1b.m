clear all
%m equally spaced points over [0,1]
m = 4; n=4;

% Vandermonde matrix t
t = zeros(m,n);
for i = 1:n
    for j = 1:m
        t(j,i) = ((j-1)/(m-1))^(n-i);
    end
end

%fliping the vandermonde matrix t to form A
A = fliplr(t)
[v_h,r_h] = house(A)

[v,r] = householder(A)

function x = backsub(R,Q,b,n)
    x = zeros(n,1);
    c = Q'*b;
    x(n) = c(n)/R(n,n);
    for i = n-1:-1:1
        h=0;
        for j=(i+1):n
            h=h+R(i,j)
        end
    end
end
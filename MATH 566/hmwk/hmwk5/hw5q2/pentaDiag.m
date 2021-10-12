function x = pentaDiag(a, b, c, d, e, f)

n = size(a,1);

x = zeros(n,1);

for i = 1:n-2
    
    li = b(i)/a(i);
    si = d(i)/a(i);
    
    a(i+1) = a(i+1) - li*c(i);
    c(i+1) = c(i+1) - li*e(i);
    b(i+1) = b(i+1) - si*c(i);
    a(i+2) = a(i+2) - si*e(i);
    
    f(i+1) = f(i+1) - li*f(i);
    f(i+2) = f(i+2) - si*f(i);
    
end

ln  = b(n-1)/a(n-1);
a(n) = a(n) - ln*c(n-1);
f(n) = f(n) - ln*f(n-1);

x(n) = f(n)/a(n);
x(n-1) = (f(n-1) - c(n-1)*x(n))/a(n-1);

for i = n-2:-1:1
    
    x(i) = (f(i) - c(i)*x(i+1) - e(i)*x(i+2))/a(i);
    
end
end

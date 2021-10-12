
% Function isolves a pentadiagonal linear system
% input : a,b,c,d,e,and f
%output : solution of the system x

function [x] = pentadiagonal(a,b,c,d,e,f)
    
    n = length(a);
    bp = zeros(n-2,1);
    dp = zeros(n-2,1);
    x  = zeros(n,1);
    for i = 1:n-2
        bp(i) = b(i)/a(i);
        dp(i) = d(i)/a(i);
        %coefficients
        a(i+1) = a(i+1) - bp(i)*c(i);
        c(i+1) = c(i+1) - bp(i)*e(i);
        b(i+1) = b(i+1) - dp(i)*c(i);
        a(i+2) = a(i+2) - dp(i)*e(i);
        %left hand side
        f(i+1) = f(i+1) - bp(i)*f(i);
        f(i+2) = f(i+2) - dp(i)*f(i);
    end
    
    a(n) = a(n) - (b(n-1)/a(n-1))*c(n-1);
    f(n) = f(n) - (b(n-1)/a(n-1))*f(n-1);
    %backward subsititution
    x(n) = f(n)/a(n);
    x(n-1) = (f(n-1) - c(n-1)*x(n))/a(n-1);
    
    for i = n-2:-1:1
        x(i) = (f(i) - c(i)*x(i+1) - e(i)*x(i+2))/a(i);
    end
end
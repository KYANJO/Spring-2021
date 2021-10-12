clear all
%close all

clc

for m =4:7
    n = 2^m;
    a = rand(n,n); b = rand(n,n);
    c_s = strass(a,b);
    C_re = a*b;
    diff =  norm(c_s) - norm(C_re)

end
fprintf('There is almost no difference in the two methods, due to very small error')

function c = strass(a,b)
nmin = 16;
[n,n] = size(a);
if n <= nmin;
c = a*b;
else
m = n/2; u = 1:m; v = m+1:n;
p1 = strass(a(u,u)+a(v,v),b(u,u)+b(v,v));
p2 = strass(a(v,u)+a(v,v),b(u,u));
p3 = strass(a(u,u),b(u,v)-b(v,v));
p4 = strass(a(v,v),b(v,u)-b(u,u));
p5 = strass(a(u,u)+a(u,v),b(v,v));
p6 = strass(a(v,u)-a(u,u),b(u,u)+b(u,v));
p7 = strass(a(u,v)-a(v,v),b(v,u)+b(v,v));
c = [p1+p4-p5+p7,p3+p5; p2+p4, p1-p2+p3+p6];
end
end
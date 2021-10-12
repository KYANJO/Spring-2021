
clear all
clc

n = 10;

J2 = 1:n-1;
J3 = 1:n-2;

a = (1:n)';

b = (-(J2+1)./3)';

c = b;

d = (-(J3+2)./6)';
e = d;

f = [1/2; 1/6; zeros(n-4,1); 1/6; 1/2];

x = pentaDiag(a, b, c, d, e, f)

A = diag(d,-2) + diag(b,-1) + diag(a) +...
    diag(c,1) + diag(e,2);

x1=A\f

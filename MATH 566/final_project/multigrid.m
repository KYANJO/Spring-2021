clear all;
close all;

% Script for testing fd2poisson over the square [a,b]x[a,b]
a = 0; b = 1;
% Laplacian(u) = f
f = @(x,y) 10*pi^2*(1+cos(4*pi*(x+2*y))-2*sin(2*pi*(x+2*y))).*exp(sin(2*pi*(x+2*y)));
% u = g on Boundary
g = @(x,y) exp(sin(2*pi*(x+2*y)));
% Exact solution is g.
uexact = @(x,y) g(x,y);
for ii =1:3
for k = 4:6
    m = 2^(k-1);

    [umg,x,y] = fd2poissonmg(f,g,a,b,m);

end
end
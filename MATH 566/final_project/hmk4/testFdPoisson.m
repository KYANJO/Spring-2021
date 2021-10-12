clear all;
close all;

% Script for testing fd2poisson over the square [a,b]x[a,b]
a = 0; b = 1;

% Laplacian(u) = f
f = @(x,y) 10*pi^2*(1+cos(4*pi*(x+2*y))-2*sin(2*pi*(x+2*y))).*exp(sin(2*pi*(x+2*y)));  
% u = g on Boundary
g = @(x,y) exp(sin(2*pi*(x+2*y)));            

%Exact solution is g.
uexact = @(x,y) g(x,y);

%Choice of the smoother
smoother = 'dampedJacobi';   %redblack for red-black-Gaus smoother
                         %dampedJacobi for Damped Jacobi smoother

for k=4:6
    m = 2^k-1; 
    h = (b-a)/(m+1);

    [u,x,y] = fd2poisson(f,g,a,b,m);

    [udst,x,y] = fd2poissondst(f,g,a,b,m);

    [umg,x,y] = fd2poissonmg(f,g,a,b,m,smoother);

end



% Plot solution
figure, set(gcf,'DefaultAxesFontSize',10,'PaperPosition', [0 0 3.5 3.5]), 
surf(x,y,udst), xlabel('x'), ylabel('y'), zlabel('u(x,y)'),
title(strcat('Numerical Solution,udst, to Poisson Equation, h=',num2str(h)));

% Plot error
figure, set(gcf,'DefaultAxesFontSize',10,'PaperPosition', [0 0 3.5 3.5]), 
surf(x,y,udst-uexact(x,y)),xlabel('x'),ylabel('y'), zlabel('Error'), 
title(strcat('Errordst, h=',num2str(h)));



% Plot solution
figure, set(gcf,'DefaultAxesFontSize',10,'PaperPosition', [0 0 3.5 3.5]), 
surf(x,y,umg), xlabel('x'), ylabel('y'), zlabel('u(x,y)'),
title(strcat('Numerical Solution,umg, to Poisson Equation, h=',num2str(h)));

% Plot error
figure, set(gcf,'DefaultAxesFontSize',10,'PaperPosition', [0 0 3.5 3.5]), 
surf(x,y,umg-uexact(x,y)),xlabel('x'),ylabel('y'), zlabel('Error'), 
title(strcat('Errormg, h=',num2str(h)));

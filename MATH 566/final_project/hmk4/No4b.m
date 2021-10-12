% Uses SOR function to to solve the poisson eqaution from problem 2 for
% various values of m and produce plots and tables that clearly show the
% forth order accuracy of the method.

a=0; b=1;

% Laplacian(u) = f
f = @(x,y) 10*pi^2*(1+cos(4*pi*(x+2*y))-2*sin(2*pi*(x+2*y))).*exp(sin(2*pi*(x+2*y)));  
% u = g on Boundary
g = @(x,y) exp(sin(2*pi*(x+2*y)));           

%Table showing the forth order acuracy of the method.
k1 = zeros(4,1);
h1=zeros(4,1);
L2=zeros(4,1);
m1=zeros(4,1);

for k = 4:7
    k1(k-3) = k;
    m1(k-3) = (2^k) - 1;
    m = (2^k) - 1;
    h1(k-3) = (b-a)/(m+1);
    h = (b-a)/(m+1);
    
    w = 2/(1+sin(pi*h)); %optimal relaxation parameter
    
    [x,y] = meshgrid(a:h:b);
    
    %Numerical solution
    [u,x,y] = SOR(f,g,a,b,m,w);
    
    % Exact solution is g.
    uexact = @(x,y) g(x,y);
   
    %Error
    error = u -uexact(x,y);
    
    %Relative 2-norm
    L2(k-3) = R2Norm(error,uexact(x,y));
    
    % Plot solution
    figure, set(gcf,'DefaultAxesFontSize',10,'PaperPosition', [0 0 3.5 3.5]),
    surf(x,y,u), xlabel('x'), ylabel('y'), zlabel('u(x,y)'),
    title(strcat('Numerical Solution to Poisson Equation, h=',num2str(h)));

    % Plot error
    figure, set(gcf,'DefaultAxesFontSize',10,'PaperPosition', [0 0 3.5 3.5]),
    surf(x,y,u-uexact(x,y)),xlabel('x'),ylabel('y'), zlabel('Error'),
    title(strcat('Error, h=',num2str(h)));

end 

%table
T = table(k1(:),m1(:),h1(:),L2(:), 'VariableNames',{'k','m','h','R2-norm'})

%polyfit
p=polyfit(log(circshift(h1,size(h1))),log(L2),1);
p
fprintf('Since the order of convergence,p, is 4.1172, which is approximately 4, \n hence the method is fourth order accurate.\n')

plot(h1,L2);
xlabel('h');
ylabel('R 2-norm');
title('A graph of h against R 2-norm');

function L2 = R2Norm(error, uexact)
    R = error .^2;
    u_ex = uexact.^2;
    L2 = sqrt(sum(R,'all')/sum(u_ex,'all'));
end
    



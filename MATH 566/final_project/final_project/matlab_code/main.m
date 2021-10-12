
%==========================================================================
% Author        : Brian Kyanjo
% Supervised by : Prof. Grady wright 
% Class         : MATH566
% Date          : March 18th, 2021
%
% Main code for calling all necessary functions. And the sample fuctions to
% be used are defined here.
% Multigrid and DST solvers are used to solve the problem in question.
% The wall time and mean wall time for each method in 3 runs has been
% performed.
%==========================================================================

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

% Compute and time the solution
k1 = zeros(1,4); h1 = zeros(1,4); m1 = zeros(1,4); 
t = zeros(1,4); 
t_dst = zeros(1,4); 
t_mg = zeros(1,4);

t1 = [];
tdst = [];
tmg = [];
for ii=1:3
    for k=4:7
    %k=10;
        
        k1(k-3) = k;
        m1(k-3) = 2^k-1;
        m = 2^k-1; 
        h1(k-3) = (b-a)/(m+1);
        h = (b-a)/(m+1);
        
        tic
        [u,x,y] = sge_solver(f,g,a,b,m);
        gedirect = toc;
        t(k-3) = gedirect;
        
        tic
        [udst,x,y] = dst_solver(f,g,a,b,m);
        gedirect = toc;
        t_dst(k-3) = gedirect;
        
        tic
        [umg,x,y] = multigrid_solver(f,g,a,b,m);
        gedirect = toc;
        t_mg(k-3) = gedirect;
    end 
t1 = [t1,t];
tdst = [tdst, t_dst]; 
tmg = [tmg,t_mg];
end

%k=4
c4=[t1(1);t1(5);t1(9)]'; 
fd4=[tdst(1);tdst(5);tdst(9)]'; 
h4=[tmg(1);tmg(5);tmg(9)]';

%k=5
c5=[t1(2);t1(6);t1(10)]'; 
fd5=[tdst(2);tdst(6);tdst(10)]'; 
h5=[tmg(2);tmg(6);tmg(10)]';

%k=6
c6=[t1(3);t1(7);t1(11)]'; 
fd6=[tdst(3);tdst(7);tdst(11)]'; 
h6=[tmg(3);tmg(7);tmg(11)]';

%k=7
c7=[t1(4);t1(8);t1(12)]'; 
fd7=[tdst(4);tdst(8);tdst(12)]'; 
h7=[tmg(4);tmg(8);tmg(12)]';

k4 = [k1(1);k1(1);k1(1)];
m4 = [m1(1);m1(1);m1(1)];
h_4 = [h1(1);h1(1);h1(1)];

%Table showing timing results of each method and for each value of m.
Table4 = table(k4,m4,h_4,c4(:),fd4(:),h4(:), 'VariableNames',{'k','m','h','t_SGE','time_DST','time_MG'});

k5 = [k1(2);k1(2);k1(2)];
m5 = [m1(2);m1(2);m1(2)];
h_5 = [h1(2);h1(2);h1(2)];
%Table showing timing results of each method and for each value of m.
Table5 = table(k5,m5,h_5,c5(:),fd5(:),h5(:), 'VariableNames',{'k','m','h','t_SGE','time_DST','time_MG'});

k6 = [k1(3);k1(3);k1(3)];
m6 = [m1(3);m1(3);m1(3)];
h_6 = [h1(3);h1(3);h1(3)];
%Table showing timing results of each method and for each value of m.
Table6 = table(k6,m6,h_6,c6(:),fd6(:),h6(:), 'VariableNames',{'k','m','h','t_SGE','time_DST','time_MG'});

k7 = [k1(4);k1(4);k1(4)];
m7 = [m1(4);m1(4);m1(4)];
h_7 = [h1(4);h1(4);h1(4)];
%Table showing timing results of each method and for each value of m.
Table7 = table(k7,m7,h_7,c7(:),fd7(:),h7(:), 'VariableNames',{'k','m','h','t_SGE','time_DST','time_MG'});


%General Table
Table = [Table4; Table5; Table6; Table7]


%mean
Tablem4 = table(k1(1),m1(1),h1(1),mean(c4),mean(fd4),mean(h4), 'VariableNames',{'k','m','h','t_SGE','time_DST','time_MG'});
Tablem5 = table(k1(2),m1(2),h1(2),mean(c5),mean(fd5),mean(h5), 'VariableNames',{'k','m','h','t_SGE','time_DST','time_MG'});
Tablem6 = table(k1(3),m1(3),h1(3),mean(c6),mean(fd6),mean(h6), 'VariableNames',{'k','m','h','t_SGE','time_DST','time_MG'});
Tablem7 = table(k1(4),m1(4),h1(4),mean(c7),mean(fd7),mean(h7), 'VariableNames',{'k','m','h','t_SGE','time_DST','time_MG'});

Table_mean = [Tablem4; Tablem5; Tablem6; Tablem7]

% 
% % Plot solution
% figure, set(gcf,'DefaultAxesFontSize',10,'PaperPosition', [0 0 3.5 3.5]), 
% surf(x,y,udst), xlabel('x'), ylabel('y'), zlabel('u(x,y)'),
% title(strcat('Numerical Solution,udst, to Poisson Equation, h=',num2str(h)));
% 
% % Plot error
% figure, set(gcf,'DefaultAxesFontSize',10,'PaperPosition', [0 0 3.5 3.5]), 
% surf(x,y,udst-uexact(x,y)),xlabel('x'),ylabel('y'), zlabel('Error'), 
% title(strcat('Errordst, h=',num2str(h)));


% %Plot solution
% figure, set(gcf,'DefaultAxesFontSize',10,'PaperPosition', [0 0 3.5 3.5]), 
% surf(x,y,umg), xlabel('x'), ylabel('y'), zlabel('u(x,y)'),
% title(strcat('Numerical Solution,umg, to Poisson Equation, h=',num2str(h)));
% 
% %Plot error
% figure, set(gcf,'DefaultAxesFontSize',10,'PaperPosition', [0 0 3.5 3.5]), 
% surf(x,y,umg-uexact(x,y)),xlabel('x'),ylabel('y'), zlabel('Error'), 
% title(strcat('Errormg, h=',num2str(h)));

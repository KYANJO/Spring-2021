clear all
close all
%clc
global IMAX 
% ADER scheme for the linear scalar advection equation 
%  q_t + a q_x = 0
% physical domain and physical parameters 
xL = -1;                % left boundary 
xR = +1;                % right boundary 
time = 0;               % initial time 
tend = 2.0;             % final time 
a = 1;                  % advection speed 
% discretization parameters 
IMAX = 200;             % number of finite volumes 
dx = (xR-xL)/IMAX;      % mesh spacing 
CFL = 0.9;              % Courant-Friedrichs-Lewy number 
NMAX = 10000;           % maximum number of time steps 
nGP  = 3;               % number of Gaussian quadrature points 
[xiGP,wGP]=gauleg(nGP);  % 
% computational mesh and initial condition 
for i=1:IMAX
    x(i) = xL + dx/2 + (i-1)*dx;   % coordinates of the cell barycenters 
    q(i) = 0;                      % compute the initial cell averages properly via numerical quadrature 
    for iGP=1:nGP
        xGP  = x(i)-dx/2 + dx*xiGP(iGP);  
    	q(i) = q(i) + wGP(iGP)*h(xGP);  
    end
end
qav0 = q;                           % save the initial cell averages, which are also the exact solution of the problem 
plot(x,q,'o')
% --------- 
% time loop 
% --------- 
tic
for n=1:NMAX
    dt = CFL*dx/abs(a);     % time step according to the CFL condition 
    if( time+dt>tend )
        dt=tend-time;       % last time step to reach the final time 
    end
    if( time>=tend )
        break               % tend reached? then stop. 
    end
    for i=1:IMAX        
        % right ADER flux 
        im1 = idx(i-1); 
        im2 = idx(i-2); 
        im3 = idx(i-3); 
        i0  = idx(i); 
        ip1 = idx(i+1); 
        ip2 = idx(i+2); 
        q0     = -(1/30)*q(ip2)-(1/60)*q(im3)+(7/60)*q(im2)-(23/60)*q(im1)+(19/20)*q(i0)+(11/30)*q(ip1); 
        qx     = -(1/180)*(2*q(im3)-10*q(im2)+5*q(im1)+205*q(i0)-215*q(ip1)+13*q(ip2))/dx;  
        qxx    = (1/8)*(q(im3)-7*q(im2)+22*q(im1)-26*q(i0)+9*q(ip1)+q(ip2))/dx^2; 
        qxxx   = (1/6)*(q(im3)-5*q(im2)+4*q(im1)+8*q(i0)-13*q(ip1)+5*q(ip2))/dx^3; 
        qxxxx  = -(1/2)*(q(im3)-7*q(im2)+18*q(im1)-22*q(i0)+13*q(ip1)-3*q(ip2))/dx^4; 
        qxxxxx = -(q(im3)-5*q(im2)+10*q(im1)-10*q(i0)+5*q(ip1)-q(ip2))/dx^5; 
        qp     = q0; 
        qt     = -a^1*qx; 
        qtt    =  a^2*qxx; 
        qttt   = -a^3*qxxx; 
        qtttt  =  a^4*qxxxx; 
        qttttt = -a^5*qxxxxx; 
        % time average of the flux of the solution of the GRP 
        fp = a*( qp*dt + qt*dt^2/2 + qtt*dt^3/6 + qttt*dt^4/24 + qtttt*dt^5/120 + qttttt*dt^6/720 )/dt; 
        % left ADER flux 
        im1 = idx(i-1-1); 
        im2 = idx(i-2-1); 
        im3 = idx(i-3-1); 
        i0  = idx(i-1); 
        ip1 = idx(i+1-1); 
        ip2 = idx(i+2-1); 
        q0     = -(1/30)*q(ip2)-(1/60)*q(im3)+(7/60)*q(im2)-(23/60)*q(im1)+(19/20)*q(i0)+(11/30)*q(ip1); 
        qx     = -(1/180)*(2*q(im3)-10*q(im2)+5*q(im1)+205*q(i0)-215*q(ip1)+13*q(ip2))/dx;  
        qxx    = (1/8)*(q(im3)-7*q(im2)+22*q(im1)-26*q(i0)+9*q(ip1)+q(ip2))/dx^2; 
        qxxx   = (1/6)*(q(im3)-5*q(im2)+4*q(im1)+8*q(i0)-13*q(ip1)+5*q(ip2))/dx^3; 
        qxxxx  = -(1/2)*(q(im3)-7*q(im2)+18*q(im1)-22*q(i0)+13*q(ip1)-3*q(ip2))/dx^4; 
        qxxxxx = -(q(im3)-5*q(im2)+10*q(im1)-10*q(i0)+5*q(ip1)-q(ip2))/dx^5; 
        qm     = q0; 
        qt     = -a^1*qx; 
        qtt    =  a^2*qxx; 
        qttt   = -a^3*qxxx; 
        qtttt  =  a^4*qxxxx; 
        qttttt = -a^5*qxxxxx; 
        fm = a*( qm*dt + qt*dt^2/2 + qtt*dt^3/6 + qttt*dt^4/24 + qtttt*dt^5/120 + qttttt*dt^6/720 )/dt; 
        qnew(i) = q(i) - dt/dx*( fp - fm ); 
    end
    time = time + dt;   % update time 
    q = qnew;           % overwrite solution 
    plot(x,q,'o')       % plot result  
    title(sprintf('Current time = %f',time)) 
    drawnow             % 
end
toc 
L2error = 0; 
for i=1:IMAX
    L2error = L2error + dx*(q(i)-qav0(i))^2; 
end
L2error = sqrt(L2error); 
disp(sprintf(' Mesh spacing dx = %e, L2error = %e ', dx, L2error)) 
xe = linspace(xL,xR,10*IMAX); 
for i=1:10*IMAX
    qe(i) = h(xe(i));  % exact solution of the problem 
end
hold on 
plot(xe,qe,'r-')
legend('ADER6','Exact') 





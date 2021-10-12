clear all
close all
%clc
global IMAX dx a 
% Runge-Kutta FV scheme for the linear scalar advection equation  
%  q_t + a q_x = 0
rk = 5; 
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
q = zeros(IMAX,1); 
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
    switch(rk)
        case(3) 
            % use the third order TVD RK scheme     
            q1   = q + dt*Lh6(q);  
            q2   = 3/4*q + 1/4*q1 + 1/4*dt*Lh6(q1);  
            qnew = 1/3*q + 2/3*q2 + 2/3*dt*Lh6(q2);  
        case(4) 
            % use the classical RK4 scheme     
            k1 = Lh6(q); 
            k2 = Lh6(q+dt/2*k1); 
            k3 = Lh6(q+dt/2*k2); 
            k4 = Lh6(q+dt*k3); 
            qnew = q + dt/6*(k1+2*k2+2*k3+k4);     
        case(5) 
            % use Butcher's RK5 scheme     
            k1 = Lh6(q); 
            k2 = Lh6(q+dt/4*k1);  
            k3 = Lh6(q+dt/8*k1+dt/8*k2);  
            k4 = Lh6(q-dt/2*k2+dt*k3); 
            k5 = Lh6(q+dt*3/16*k1+dt*9/16*k4); 
            k6 = Lh6(q-dt*3/7*k1+dt*2/7*k2+dt*12/7*k3-dt*12/7*k4+dt*8/7*k5); 
            qnew = q + dt/90*(7*k1+32*k3+12*k4+32*k5+7*k6);  
    end
    time = time + dt;   % update time 
    q = qnew;           % overwrite solution 
    plot(x,q,'ko')      % plot result  
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
legend('RKFV6','Exact') 





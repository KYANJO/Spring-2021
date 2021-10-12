clear all
close all
%clc
global IMAX dx a
% MUSCL-Hancock TVD FV scheme for the linear scalar advection equation
%  q_t + a q_x = 0
% physical domain and physical parameters
xL = -1;                % left boundary
xR = +1;                % right boundary
time = 0;               % initial time
tend = 10.0;            % final time
a = 1;                  % advection speed
% discretization parameters
IMAX = 100;             % number of finite volumes
dx = (xR-xL)/IMAX;      % mesh spacing
CFL = 0.9;              % Courant-Friedrichs-Lewy number
NMAX = 10000;           % maximum number of time steps
% computational mesh and initial condition
q = zeros(IMAX,1);
for i=1:IMAX
    x(i) = xL + dx/2 + (i-1)*dx;   % coordinates of the cell barycenters
    q(i) = h(x(i));                % compute the initial cell averages using the midpoint rule
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
    % slope reconstruction
    for i=1:IMAX
        im1 = idx(i-1);
        ip1 = idx(i+1);
        % calculate the limited slope using the minmod limiter
        slope = minmod( (q(ip1)-q(i))/dx, (q(i)-q(im1))/dx );
        % calculate the boundary-extrapolated values
        wR(i) = q(i) + 0.5*dx*slope;
        wL(i) = q(i) - 0.5*dx*slope;
        % time evolution 
        qt = -a*slope; 
        wR(i) = wR(i) + 0.5*dt*qt; 
        wL(i) = wL(i) + 0.5*dt*qt; 
    end
    % flux calculation using the upwind scheme
    for i=1:IMAX
        im1 = idx(i-1);
        ip1 = idx(i+1);
        fp = 0.5*a*(wR(i)+wL(ip1))-0.5*abs(a)*(wL(ip1)-wR(i));
        fm = 0.5*a*(wR(im1)+wL(i))-0.5*abs(a)*(wL(i)-wR(im1));
        qnew(i) = q(i) - dt/dx*( fp - fm );
    end
    
    
    time = time + dt;   % update time
    q = qnew;           % overwrite solution
    plot(x,q,'ko')      % plot result
    title(sprintf('Current time = %f',time))
    axis([xL xR 0 1]) 
    drawnow              
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
legend('MUSCL','Exact')





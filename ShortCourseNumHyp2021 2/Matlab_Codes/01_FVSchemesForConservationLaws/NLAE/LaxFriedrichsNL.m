clear all
close all
clc
% FV method for the solution of the PDE q_t+f_x=0 using Lax-Friedrichs flux

% Domain definition
xL = -1;
xR = 1;
time = 0; tend = 0.4;
% Number of nodes
IMAX = 100;
% Discretization space
dx = (xR-xL)/IMAX;
x = linspace(xL+dx/2,xR-dx/2,IMAX);
% Maximum number of time steps
NMAX = 100000;
% CFL
CFL = 0.9;

% Initial condition
qL = 1; qR= 2; % Initial states
q = zeros(IMAX,1); qnew = zeros(IMAX,1);
for i=1:IMAX
    if (x(i)<=0)
        q(i)=qL;
    else
        q(i) = qR;
    end
end
% Plot initial solution
plot(x,q,'r')

% Time loop 
for  n =1:NMAX
    % Compute time step
    amax = 0;
    for i=1:IMAX
        amax = max( amax, abs( a( q(i) ) ) );
    end
    dt = CFL*dx/amax;
    % Ensure exact tend
    if (time +dt > tend)
        dt = tend -time;
    end
    % Stop criteria
    if (time>=tend)
        break
    end    
    
    % FV method
    for  i = 1:IMAX
        if (i==1)
            qnew(i) = qL;
        elseif (i==IMAX)
            qnew(i) = qR;
        else
            % Lax-Friedrichs flux
            fp = 0.5*(f(q(i+1)) + f(q(i))) - 0.5*dx/dt*(q(i+1)-q(i) ); % f_{i+1/2}
            fm = 0.5*(f(q(i)) + f(q(i-1))) - 0.5*dx/dt*(q(i)-q(i-1) ); % f_{i-1/2}
            % Solution at n+1
            qnew(i) = q(i) -dt/dx*(fp-fm);
        end
    end
    % Overwrite solution
    q = qnew;
    % Update time step
    time = time+dt;
    % Plot
    plot(x,q,'ro')
    title(sprintf('Time %f',time))
    drawnow
end

% Plot exact solution
xe = linspace(xL,xR,10*IMAX);
for i=1:length(xe)
    xi = xe(i)/time;
    qe(i)= ExactRiemannSolverNL(qL,qR,xi);    
end
hold on
plot(xe,qe,'k-')
legend('Lax-Friedrichs','Exact solution')
hold off







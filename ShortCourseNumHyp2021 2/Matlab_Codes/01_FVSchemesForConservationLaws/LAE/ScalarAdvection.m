% Upwind method to solve the scalar advection equation
% q_t + a q_x = 0
clear all
close all
clc

% Physical setup
xL = -1;    % Left boundary
xR = 1;     % Right boundary
t = 0;      % Initial time
tend = 0.5; % End time
a =1;       % Transport velocity

% Define data related to the mesh
imax = 1000;            % Number of nodes
dx = (xR-xL)/(imax-1); 	% Cell length
CFL = 0.9;              % Courant-Friedrichs-Lewy number
Nmax = 10000;           % MAximum number of time steps

for i = 1:imax
    % Define mesh
    x(i) = xL + (i-1)*dx;
    % Compute IC in the mash nodes
    q(i) = h(x(i));
end
plot(x,q,'o')

% Numerical method
for  n=1:Nmax
    % Compute time step
    dt = CFL*dx/abs(a);
    % Ensure that we exactly attain the final time
    if (t+dt>tend)
        dt = tend-t;
    end
    % Stop criteria
    if (t>=tend)
        break
    end
    
    % Calcolo della soluzione al nuovo passo di tempo per ogni nodo
    for i =1:imax
        if (i==1)
            % Left boundary conditions
            qnew(i) = h(xL-a*(t+dt));
        elseif (i==imax)
            % Right boundary conditions
            qnew(i) = h(xR-a*(t+dt));
        else       
            % Upwind method
            if (a>0) 
                % Upwind method per a>0
                qnew(i) = q(i) - a*dt/dx*(q(i)-q(i-1)); % Bacward-upwind / backward-downwind
            else
                % Upwind method per a<0
                 qnew(i) = q(i) - a*dt/dx*(q(i+1)-q(i)); % forward-upwind / forward-downwind
            end
            % Centred
            %qnew(i) = q(i) - 0.5*a*dt/dx*(q(i+1)-q(i-1)); 
            % Lax-Friedricchs
            %qnew(i) = 0.5*(q(i+1)+q(i-1)) - 0.5*a*dt/dx*(q(i+1)-q(i-1)); 
            % Lax-Wendroff
            %             qnew(i) = q(i) - 0.5*a*dt/dx*(q(i+1)-q(i-1)) ...
            %                 + 0.5*(a*dt/dx)^2 * (q(i+1)-2*q(i)+q(i-1)); 
            %             % FORCE
            %             qnew(i) = 0.5*(0.5*(q(i+1)+q(i-1)) - 0.5*a*dt/dx*(q(i+1)-q(i-1)) + ...
            %                 q(i) - 0.5*a*dt/dx*(q(i+1)-q(i-1)) ...
            %                 + 0.5*(a*dt/dx)^2 * (q(i+1)-2*q(i)+q(i-1)));
        end
    end
    
    % Update time step
    t = t + dt;
    % Overwrite the solution
    q = qnew;
    % Plot
    plot(x,q,'o')
    % Force plotting now
    drawnow   
end
%% Plot exact solution
imaxe = 100000; % Number of points for the discreetization of the exact solution
xe = linspace(xL,xR,imaxe); % Griglia per la soluzione essata
for i=1:imaxe
    qe(i) = h(xe(i)-a*t);
end
hold on
plot(xe,qe,'k-')
legend('Upwind','Exact solution')
hold off
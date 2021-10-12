% Code to solve Shallow water equations with variable bottom
% using FV methods with Rusanov scheme
%**************************************************************************
close all
clear all
clc
global g

%% Define model parameters
g = 9.81;           % Gravity
Ks = 26;            % Strickler coefficient

%% Define computational domain and compute grid
xL = -1;            % Left boundary
xR = 1;             % Right boundary
IMAX = 100;         % Number of cells
dx = (xR-xL)/IMAX;  % Mesh spacing
x = linspace(xL+dx/2,xR-dx/2,IMAX); % Cell baricentres

%% Time
time = 0;           % Initial time
tend = 10;          % End time
CFL = 0.9;          %Courant number
NMAX = 1000;        % Maximum number of time steps

%% Initial condition (lake at rest)
for i = 1:IMAX
    eta = 2;        % Free surface
    sigma = 0.1;    % Gaussian halfwidth
    b = exp(-0.5*x(i)^2/sigma^2);   % Bottom bathymetry
    if (x(i)>0.5)
        b = b+0.5;
    end
    h = eta - b;    % water depth
    u = 0;          % Velocity
    Q(:,i) = [h;h*u;b];  % Conservative variables
end
% Plot IC
plot(x,Q(1,:)+Q(3,:),'bo')
hold on
plot(x,Q(3,:),'ko')

%% Main solver, time loop
for n=1:NMAX
    % Compute time step
    amax = 0;
    for i=1:IMAX
        L = lambda(Q(:,i));
        amax = max(amax,max(abs(L))); 
    end
    dt = CFL*dx/amax;
    if (time +dt >tend)
        dt = tend-time;
    end
    if (time>=tend)
        break
    end
    % Compute solution at the new time
    for i=1:IMAX
        if (i==1)
            [fp,Dp] = RusanovFluxWB(Q(:,i),Q(:,i+1));
            QBC = Q(:,i);
            QBC(2) = -QBC(2);
            [fm,Dm] = RusanovFluxWB(QBC,Q(:,i));
        elseif (i==IMAX)
            QBC = Q(:,i);
            QBC(2) = -QBC(2);
            [fp,Dp] = RusanovFluxWB(Q(:,i),QBC);
            [fm,Dm] = RusanovFluxWB(Q(:,i-1),Q(:,i));            
        else
            [fp,Dp] = RusanovFluxWB(Q(:,i),Q(:,i+1));
            [fm,Dm] = RusanovFluxWB(Q(:,i-1),Q(:,i));
        end
        Qnew(:,i) = Q(:,i) -dt/dx*(fp-fm) -dt/dx*(Dp+Dm);
    end
    % Update time and overwrite solution
    Q = Qnew;
    time = time + dt;    
    % Plot
    subplot(2,1,1)
    hold off
    plot(x,Q(1,:)+Q(3,:),'bo')
    ylabel('eta, b')
    hold on
    plot(x,Q(3,:),'ko')
    title(sprintf('Current time = %f',time))
    subplot(2,1,2)
    plot(x,Q(2,:)./Q(1,:))
    ylabel('u')
    drawnow
end





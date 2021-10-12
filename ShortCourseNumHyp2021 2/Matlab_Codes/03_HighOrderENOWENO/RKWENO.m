%% Matlab code for the solution of  non-linear hyperbolic systems
% Q_t + F_x(Q) + B(Q)Q_x = 0  (PDE)
% Using Ruge-Kutta WENO finite volume schemes
clear all
close all
clc

global gamma
global N IMAX nVar dx dt 
global pdetype flux

%% Select the model
% 0 == Euler equations
pdetype = 0;
gamma = 1.4;

%% Select method
% Flux
% 0 == Rusanov flux
flux = 0;
% Order of the final WENO scheme
% 2 == TVD schme with minmod
% 3 == WENO3
% 5 == WENO5
N = 3;
% Runge-Kutta type
rk= 5
% CFL
CFL = 0.9;

% Number of variables
switch(pdetype)
    case(0)
        nVar = 3;
end

%% Define the computational domain and grid
x0 = -1;    % Left boundary
x1 = 1;     % right boundary
IMAX = 200; % Number of cells
dx = (x1-x0)/IMAX; 
x = linspace(x0+dx/2,x1-dx/2,IMAX);
time = 0;   % Initial time
tend = 0.4; % End time
NMAX = 100000; % Maximum number of time iterations
nGP = 3;    % Number of Gauss points
[xiGP,wGP] = gauleg(nGP); % Coordinates and weights for Gauss-Legendre quadrature points

%% Define initial solution (primitive variables)
switch(pdetype)
    case(0) % Sod test
        VL = [1;0;1];           % Left state
        VR = [0.125; 0; 0.1];   % Right state
        tend = 0.2;             % Final time
        xc = 0;                 % Position of the initial discontinuity
end
uL = PDEPrim2Cons(VL);
uR = PDEPrim2Cons(VR);

%% Project the IC into the spacial basis functions
q = zeros(nVar,IMAX);
for i=1:IMAX
    for iGP = 1:nGP
        %j = j+1;
        xi = xiGP(iGP);
        xGP = x(i) - dx/2 + xi*dx;
        % Initialize assuming that we have a Riemann problem as IC
        if (xGP<xc)
            u0 = uL;
        else
            u0 = uR;
        end
        weight = wGP(iGP); % Weight of the G point
        q(:,i) = q(:,i) + weight*u0(:);
    end
end

% Recover the names of the variables
varname = PDEVarName;
qnew = q;

%% Time loop
disp('Preparation done. Starting code...')
for n=1:NMAX
    % Compute time step
    lmax = 0;
    for i=1:IMAX
        lambda = max( abs( PDEEigenvalues(q(:,i)) ) );
        lmax = max(lmax,lambda);
    end
    dt = CFL*dx/lmax;
    % Ensure arriving to the final time
    if (time +dt > tend)
        dt = tend-time;
    end
    % Stop criterium
    if (time >= tend)
        break
    end
    
    % Runge-Kutta
    switch(rk)
        case(3)
            % Third order TVD RK 
            q1 = q + dt*LhWENO(q); 
            q2 = 3/4*q + 1/4*q1 + 1/4*dt*LhWENO(q1);
            qnew = 1/3*q + 2/3*q2 + 2/3*dt*LhWENO(q2);
        case(4)
            % Classical RK4 scheme
            k1 = LhWENO(q);
            k2 = LhWENO(q+dt/2*k1);
            k3 = LhWENO(q+dt/2*k2);
            k4 = LhWENO(q+dt*k3);
            qnew = q +dt/6*(k1+2*k2+2*k3+k4);
        case(5) 
            % Butcher's RK5 scheme     
            k1 = LhWENO(q); 
            k2 = LhWENO(q+dt/4*k1);  
            k3 = LhWENO(q+dt/8*k1+dt/8*k2);  
            k4 = LhWENO(q-dt/2*k2+dt*k3); 
            k5 = LhWENO(q+dt*3/16*k1+dt*9/16*k4); 
            k6 = LhWENO(q-dt*3/7*k1+dt*2/7*k2+dt*12/7*k3-dt*12/7*k4+dt*8/7*k5); 
            qnew = q + dt/90*(7*k1+32*k3+12*k4+32*k5+7*k6); 
    end
    % Update time and overwrite solutio
    q = qnew;
    time = time+dt;
    for i=1:IMAX
        v(:,i) = PDECons2Prim(q(:,i));
    end
    
    % Plot the solution
    for iVar=1:nVar
        subplot(nVar,1,iVar)
        plot(x,v(iVar,:),'o');
        xlabel('x')
        ylabel(varname(iVar))
    end
    title(sprintf('Time = %f',time))
    drawnow
end

%% Plot exact solution
switch(pdetype)
    case(0) % Euler
        xe = linspace(x0,x1,10*IMAX);
        for i=1:length(xe)
            s = (xe(i)-xc)/time;
            [rho,uu,pp] = ExactRiemannEuler(VL(1),VR(1),VL(2),VR(2),VL(3),VR(3),s,gamma,gamma);
            ue(1,i) = rho;
            ue(2,i) = uu;
            ue(3,i) = pp;
        end
        for iVar=1:nVar
            subplot(nVar,1,iVar);
            hold on
            plot(xe,ue(iVar,:),'r-')
        end
end



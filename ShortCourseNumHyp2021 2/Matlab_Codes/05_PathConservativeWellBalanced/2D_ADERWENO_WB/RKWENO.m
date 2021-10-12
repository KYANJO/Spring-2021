%% Matlab code for the solution of the non-linear hyperbolic system
%     Q_t +f_x(Q) +B(q) Q_x = 0 (PDE)
% using a Runge-Kutta WENO finite volume scheme
%*******************************************************************************
clear all
close all
clc

global pdetype flux rk
global gamma g
global N IMAX nVar dx dt IWB

%% Select model to be solved
% 0 == Euler equations
% 1 == SW with scalar transport
pdetype = 1;
% Equation parameters
gamma = 1.4;            % Ratio of specific heats for Euler
g           = 9.81;          % Gravity constant for SW 

%% Select method
% Flux: 
% 0 == Rusanov; 
% 1 == Osher-type (DOT)
flux = 0;
% Order of the final WENO scheme 
%  2 == TVD scheme with minmod 
%  3 == WENO3
%  5 == WENO5 
N = 2;
% Runge-Kutta type
rk = 3;
% Courant number
CFL = 0.9;

%% Define grid, time and Gauss points
x0         = -1;                        % left boundary location
x1         = 1;                          % Right boundary location
IMAX    = 200;                      % number of elements
dx = (x1-x0)/IMAX;              % space increment
x = linspace(x0+dx/2,x1-dx/2,IMAX);     % nodes
time = 0;                               % Initial time
tEnd = 0.4;                            % final time
NMAX  = 10000;                  % max number of time iterations
nGP  = 3;                              % Number of Gauss points
[xiGP,wGP]=gauleg(nGP);  % Coordinates and weights

%% Define initial solution in primitive variables
switch(pdetype)
    case(0)
        % Sod problem
        VL = [1;0;1];               % Left state
        VR = [0.125,0,0.1];     % Right state
        tEnd = 0.2;
        xc = 0;                        % Initial position of the discontinuity
    case(1)
        VL = [2;0;1];  
        VR = [1;0;0];
        tEnd  = 0.15;
         xc = 0;                        % Initial position of the discontinuity
end
% Compute IC in conservative variables
uL = PDEPrim2Cons(VL);
uR = PDEPrim2Cons(VR);

%% Recall number of variables attending to the model
switch(pdetype)
    case(0)
        nVar  = 3;
    case(1)
        nVar  = 3;
end

%% Project IC onto the spatial basis functions
disp('Projecting initial condition')
q = zeros(nVar,IMAX);
for i=1:IMAX
    for ixGP = 1:nGP
        j = j +1;
        xi = xiGP(ixGP);
        xGP = x(i) -dx/2 + xi*dx;
        % Initialize assuming we have a Riemann problem as IC
        % Change here if necessary
        if (xGP<xc)
            u0 = uL;
        else
            u0 = uR;
        end
        weight = wGP(ixGP);
        q(:,i) = q(:,i) + weight*u0(:);
    end
end

%% Main time loop
disp('Preparation done. Starting code...')
for n = 1:NMAX
    %% Compute time step
    lmax = 0;
    for i=1:IMAX
        lambda = max(abs(PDEEigenvalues(q(:,i))));
        lmax = max(lmax,lambda);
    end
    dt = CFL*dx/lmax;
    % Stop criterium
    if (time>=tEnd)
        break
    end
    % Final time step
    if (time+dt>tEnd)
        dt = tEnd-time;
    end
    
    %% Runge-Kutta
    switch(rk)
        case(3)
            % Third order TVD RK scheme 
            q1 = q + dt*LhWENO(q);
            q2   = 3/4*q + 1/4*q1 + 1/4*dt*LhWENO(q1);  
            qnew = 1/3*q + 2/3*q2 + 2/3*dt*LhWENO(q2);  
        case(4)
            % Classical RK4 scheme
            k1 = LhWENO(q); 
            k2 = LhWENO(q+dt/2*k1); 
            k3 = LhWENO(q+dt/2*k2); 
            k4 = LhWENO(q+dt*k3); 
            qnew = q + dt/6*(k1+2*k2+2*k3+k4);
        case(5)
            % use Butcher's RK5 scheme     
            k1 = LhWENO(q); 
            k2 = LhWENO(q+dt/4*k1);  
            k3 = LhWENO(q+dt/8*k1+dt/8*k2);  
            k4 = LhWENO(q-dt/2*k2+dt*k3); 
            k5 = LhWENO(q+dt*3/16*k1+dt*9/16*k4); 
            k6 = LhWENO(q-dt*3/7*k1+dt*2/7*k2+dt*12/7*k3-dt*12/7*k4+dt*8/7*k5); 
            qnew = q + dt/90*(7*k1+32*k3+12*k4+32*k5+7*k6);  
    end
    
    %% Update time and overwrite solution
    time = time +dt;
    q = qnew;
    for i=1:IMAX
        v(:,i) = PDECons2Prim(q(:,i));
    end
    
    %% Plot solution
    varname = PDEVarName;
    for iVar =1:nVar
        subplot(nVar,1,iVar);
        plot(x,v(iVar,:),'o');
        xlabel('x')
        ylabel(varname(iVar,:))
    end
    title(sprintf('Regular output at time t=%f',time))
    drawnow
end
% End of main loop

%% Plot exact solution if available
switch(pdetype)
    case(0)
        xe = linspace(x0,x1,2000);
        for i = 1:length(xe)
            s = (xe(i)-xc)/time;
            [rho,uu,pp] = ExactRiemannEuler(VL(1),VR(1),VL(2),VR(2),VL(3),VR(3),s,gamma,gamma);
            ue(1,i) = rho;
            ue(2,i) = uu;
            ue(3,i) = pp;
        end
        for iVar = 1:nVar
            subplot(nVar,1,iVar);
            hold on
            plot(xe,ue(iVar,:),'r-');
        end
    case(1)
        xe = linspace(x0,x1,2000);
        for i = 1:length(xe)
            s = (xe(i)-xc)/time;
            [dd,uu,psi] = ExactRiemannShallowWater(uL(1),uR(1),uL(2)/uL(1),uR(2)/uR(1),uL(3)/uL(1),uR(3)/uR(1),g,s);
            ue(1,i) = dd;
            ue(2,i) = uu;
            ue(3,i) = psi;
        end
        for iVar = 1:nVar
            subplot(nVar,1,iVar);
            hold on
            plot(xe,ue(iVar,:),'r-');
        end
end















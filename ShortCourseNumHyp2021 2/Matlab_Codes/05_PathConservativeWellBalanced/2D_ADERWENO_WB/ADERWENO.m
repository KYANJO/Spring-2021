%% Matlab code for the solution of the non-linear hyperbolic system
%     Q_t +f_x(Q) +B(q) Q_x = 0 (PDE)
% using an ADER WENO finite volume scheme with the new local
% space-time predictor instead of the classical Cauchy-Kovalevskaya
% procedure.
%*******************************************************************************
% For details see: S. Busto, S. Chiocchetti, M. Dumbser, E. Gaburro, I. Peshkov. 
%  High order ADER schemes for continuum mechanics. 
%  Frontiers in Physics, 8:32, 2020. DOI: 10.3389/fphy.2020.00032 
%*******************************************************************************

clear all
%close all
clc

global pdetype flux CWENO
global gamma g
global N IMAX nVar

%% Select model to be solved
% 0 == Euler equations
% 1 == SW with scalar transport
% 2 == SW with variable bottom
pdetype = 2;
% Equation parameters
gamma = 1.4;            % Ratio of specific heats for Euler
g           = 9.81;          % Gravity constant for SW 

%% Select method
% Flux: 
% 0 == Rusanov; 
% 1 == Osher-type (DOT)
flux = 0;
% Reconstruction: 
% 0 == polynomial ENO-type WENO
% 1 == polynomial CWENO
CWENO = 1;
% Polynomial degree (N=1,..,5)
N = 3;
% Courant number
CFL = 0.9;

%% Define grid, time and DOFS
x0         = -1;                        % left boundary location
x1         = 1;                          % Right boundary location
IMAX    = 200;                      % number of elements
dx = (x1-x0)/IMAX;              % space increment
x = linspace(x0+dx/2,x1-dx/2,IMAX);     % nodes
time = 0;                               % Initial time
tEnd = 0.4;                            % final time
NMAX  = 10000;                  % max number of time iterations
NFINE   = 10 + 2*(N+1);     % fine output of the results using the DG sub-cell resolution
nDOFs  = N+1;                    % number of spatial degrees of freedom
nDOF    = (N+1)^2;             % number of space-time degrees of freedom

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
    case(2)
        VL = [2;0;0]; 
        VR = [1;0;0.2];  
        tEnd  = 0.15;
         xc = 0;                        % Initial position of the discontinuity
end
% Compute IC in conservative variables
uL = PDEPrim2Cons(VL);
uR = PDEPrim2Cons(VR);

%% Initialization related to the predictor 
Initialization;

%% Build matrix for well balance attending to the system
switch(pdetype)
    case(0) % Euler
        nVar  = 3;
        IWB = eye(nVar,nVar);      
    case(1)
        nVar  = 3;
        IWB = eye(nVar,nVar); 
    case(2)
        nVar  = 3;
        % Well-balanced identity matrix 
        IWB = eye(nVar,nVar);
        IWB(1,3)=1; 
        IWB(3,3)=0;  
end

%% Project IC onto the spatial basis functions
disp('Projecting initial condition')
uhat = zeros(nVar,IMAX);
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
        uhat(:,i) = uhat(:,i) + weight*u0(:);
    end
end

%% Main time loop
disp('Preparation done. Starting code...')
for n = 1:NMAX
    %% Compute time step
    lmax = 0;
    for i=1:IMAX
        lambda = max(abs(PDEEigenvalues(uhat(:,i))));
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
    %% WENO reconstruction to obtain the high order 
    % polynomial w from the cell averages
    what = Reconstruction(uhat);
    
    %% Local space-time DG predictor
    qhat = zeros(nDOF,nVar,IMAX);
    for i=1:IMAX
        % Simple initial guess equal to the cell average
        for iVar = 1:nVar
            qhat(:,iVar,i) = what(1,iVar,i);
        end
        % Iterative method
        for iter = 1:N+2
            % 1.- Compute the derivative of qh with respect to xi (dq/dx)
            dqdxi = dudx*qhat(:,:,i);
            % 2.- Compute the flux f(Q) and the non-conservative product B*Qx
            for k= 1:nDOF
                fhat(k,:,i)   = PDEFlux(qhat(k,:,i))';
                BQxhat(k,:,i) = PDEMatrixB(qhat(k,:,i))*dqdxi(k,:)'*dt/dx;   % not necessary for Euler
            end
            % 3.- Compute the flux f*(Q) in the reference space
            fhatstar = dt/dx*fhat(:,:,i);
            % 4.- Now do the ADER / local space-time DG iteration
            qhat(:,:,i) = iK1*( F0*what(:,:,i) - Kxi*fhatstar(:,:) - M*BQxhat(:,:,i) );
        end
    end
    
    %% Compute the left and right time-integrated states and fluxes
    % This time integration is done analytically beforehand, using Fp and
    % Fm (quadrature-free implementation)
    for i = 1:IMAX
        FR(:,i) = Fp*fhat(:,:,i);
        FL(:,i) = Fm*fhat(:,:,i);
        QR(:,i) = Fp*qhat(:,:,i);
        QL(:,i) = Fm*qhat(:,:,i);
    end
    
    %% Quadrature-free one-step ADER-FV scheme
    for i = 2:IMAX-1
        ip1 = i +1;
        im1 = i-1;
        if (flux == 0) % Rusanov flux
            % Use a simple Rusanov flux at the element interfaces
            % (smax is the maximum signal speed in the system)
            % !! be careful: to compute smax here we must use what(1,:,i)
            % !! instead of uhat(:,i), since uhat is immediately overwritten
            % !! using uhat would lead to a non-conservative scheme !! 
            % Right
            smax = max( max( abs( PDEEigenvalues(what(1,:,i)) ) ), max( abs( PDEEigenvalues(what(1,:,ip1)) ) ) );
            fp = 0.5*(FL(:,ip1)+FR(:,i)) - 0.5*smax*IWB*(QL(:,ip1)-QR(:,i));
            BRoe = RoeMatrix(QR(:,i),QL(:,ip1)); % Not necessary for Euler
            Dp = 0.5*BRoe*(QL(:,ip1)-QR(:,i)); % Not necessary for Euler     
            % Left
            smax = max( max( abs( PDEEigenvalues(what(1,:,im1)) ) ), max( abs( PDEEigenvalues(what(1,:,i)) ) ) );
            fm = 0.5*(FL(:,i)+FR(:,im1)) - 0.5*smax*IWB*(QL(:,i)-QR(:,im1));
            BRoe = RoeMatrix(QR(:,im1),QL(:,i)); % Not necessary for Euler
            Dm  = 0.5*BRoe*(QL(:,i)-QR(:,im1)); % Not necessary for Euler
        else % Osher-type flux
            % Use the Osher-Solomon-type Riemann solver of Dumbser and Toro
            % M. Dumbser and E.F. Toro. A Simple Extension of the Osher Riemann Solver to Non-Conservative Hyperbolic Systems. Journal of Scientific Computing, 48:70-88, 2011 
            % M. Dumbser and E.F. Toro. On Universal Osher-Type Schemes for General Nonlinear Hyperbolic Conservation Laws. Communications in Computational Physics, 10:635-671, 2011
            % Right
            [BRoe,absA] = OsherMatrix( QR(:,i),QL(:,ip1) );
            fp = 0.5*(FL(:,ip1)+FR(:,i))-0.5*absA*(QL(:,ip1)-QR(:,i));
            Dp = 0.5*BRoe*(QL(:,ip1)-QR(:,i));
            % Left
            [BRoe,absA] = OsherMatrix( QR(:,im1),QL(:,i) );
            fm = 0.5*(FL(:,i)+FR(:,im1))-0.5*absA*(QL(:,i)-QR(:,im1));
            Dm = 0.5*BRoe*(QL(:,i)-QR(:,im1));
        end
        % Final ones-step finite volume scheme
        uhat(:,i) = uhat(:,i) - dt/dx*(fp-fm) - dt/dx*(Dm+Dp) - (MSrc(1,:)*BQxhat(:,:,i))';
    end
    
    %% Update time and overwrite solution
    time = time +dt;
    for i=1:IMAX
        vhat(:,i) = PDECons2Prim(uhat(:,i));
    end
    
    %% Plot solution
    varname = PDEVarName;
    for iVar =1:nVar
        subplot(nVar,1,iVar);
        plot(x,vhat(iVar,:),'o');
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
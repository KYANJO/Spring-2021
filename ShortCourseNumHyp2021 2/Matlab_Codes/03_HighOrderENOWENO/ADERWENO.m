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
close all
clc

global gamma
global N IMAX nVar 
global pdetype flux CWENO

%% Select the model
% 0 == Euler equations
pdetype = 0;
gamma = 1.4;

%% Select method
% Flux
% 0 == Rusanov flux
flux = 0;
% Reconstruction 
% 0 == polynomial ENO-type
% 1 == CWENO
CWENO = 1;
% Polinomial degree
N = 2;
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
tend = 0.2; % End time
NMAX = 100000; % Maximum number of time iterations
nDOFs = N+1;   % Number of degrees of freedom in space
nDOF = (N+1)^2; % Numeber of degrees of freedom in space-time
NFINE = 10+2*(N+1);

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

% Initialization
Initialization;

% Define number of variables (+ well balance matrix)
switch(pdetype)
    case(0) % Euler
        nVar = 3;
end

%% Project IC onto the spacial basis functions
uhat = zeros(nVar,IMAX);
for i=1:IMAX
    for ixGP = 1:nGP
        xi = xiGP(ixGP);
        xGP = x(i)-dx/2+xi*dx;
        % Initialize assuming the we have a RP
        if (xGP<xc)
            u0 = uL;
        else
            u0 = uR;
        end
        weight = wGP(ixGP);
        uhat(:,i) =  uhat(:,i) +weight*u0(:);
    end
end

disp('Preparation done. Starting code...')
%% Main time loop
for n=1:NMAX
    % Compute time step
    lmax = 0;
    for i=1:IMAX
        lambda = max( abs( PDEEigenvalues(uhat(:,i)) ) );
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
    
    %% 1.- WENO reconstruction
    what = Reconstruction(uhat);
    
    %% 2.- Local space-time predictor
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
            % 2.- Compute flux (and non-conservative product)
            for k =1:nDOF
                fhat(k,:,i) = PDEFlux(qhat(k,:,i))';
            end
            % 3.- Compute the flux in the reference space
            fhatstar = dt/dx*fhat(:,:,i);
            % 4.- Now do the ADER/local space-time DG iteration
            qhat(:,:,i) = iK1*(F0*what(:,:,i) - Kxi*fhatstar(:,:)); %...            
        end
    end
    
    %% 3.- Compute the left and right time-integrated states and fluxes
    for i = 1:IMAX
        FR(:,i) = Fp*fhat(:,:,i);
        FL(:,i) = Fm*fhat(:,:,i);
        QR(:,i) = Fp*qhat(:,:,i);
        QL(:,i) = Fm*qhat(:,:,i);
    end
    
    
    %% 4.- Quadrature-free one-step ADER-FV scheme
    for i = 2:IMAX-1
        ip1 = i+1;
        im1 = i-1;
        if (flux ==0) % Rusanov
            % Right
            smax = max( max(abs(PDEEigenvalues(what(1,:,i)))), ...
                max(abs(PDEEigenvalues(what(1,:,ip1))))); 
            fp = 0.5*(FL(:,ip1)+FR(:,i)) - 0.5*smax*(QL(:,ip1)-QR(:,i));
            % Left
            smax = max( max(abs(PDEEigenvalues(what(1,:,im1)))), ...
                max(abs(PDEEigenvalues(what(1,:,i))))); 
            fm = 0.5*(FL(:,i)+FR(:,im1)) - 0.5*smax*(QL(:,i)-QR(:,im1));
        %else % Osher-type flux
        end
        % 4.2.- Compute the solution using FV
        uhat(:,i) = uhat(:,i) -dt/dx*(fp-fm);
    end
      
    
    %% Update the time 
    time = time +dt;
    for i=1:IMAX
        vhat(:,i) = PDECons2Prim(uhat(:,i));
    end
    
    %% Plot the solution 
    varname = PDEVarName;
    for iVar=1:nVar
        subplot(nVar,1,iVar);
        plot(x,vhat(iVar,:),'o')
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

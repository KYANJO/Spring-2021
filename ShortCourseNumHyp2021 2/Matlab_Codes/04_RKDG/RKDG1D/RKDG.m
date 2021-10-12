% Runge-Kutta DG scheme for the scalar conservation law
% q_t +f_x = 0
%*************************************************************
clear all           % empty the memory (remove existing variables) 
close all           % closes all open figures
clc                    % clear screen

global N;           % polynomial degree of the basis functions 
global IMAX;        % number of variables in the PDE 
global nGP;         % number of Gaussian quadrature points 
global wGP;         % Gaussian weights 
global xiGP;        % Gaussian points 
global Me;          % element mass matrix 
global iMe;         % inverse of the element mass matrix 
global dx;          % mesh spacing 
global phiL phiR phiGP phi_xiGP % basis functions stored on the boundaries and in the Gaussian quadrature points 

%% Define computational domain, times and compute grid
xL = -1;            % left boundary of the domain 
xR = +1;            % right boundary of the domain
N  = 4;             % polynomial approximation degree 
IMAX = 75;          % number of mesh points  
time = 0;           % Initial time
tend = 0.5;         % Final time of the simulation
NMAX = 10000;       % Maximum number of time steps

%% Select initial condition
% 0 == Riemann problem
% 1 == Smooth function h(x)
ictype = 1;
% Left and right states for Riemann problem and initial position of the
% discontinuity
qL = 1;
qR = 2;
xc = -0.5;

%% Select scheme
% 3 == third order TVD RK
% 4 == classical RK
% 5 == Butcher's fifth order scheme  
rktype = 4;

%% Define mesh, CFL, DOF ...
dx = (xR-xL)/IMAX;                  % Mesh spacing
x = linspace(xL+dx/2,xR-dx/2,IMAX); % IMAX equidistant points
nSub = 2*(N+1);                     % Subgrid for plotting
nDOF = N+1;                         % Number of degrees of freedom per element
nGP = N+1;                          % Number of Gaussian quadrature points
[xiGP,wGP] = gauleg(nGP);           % Compute Gaussian points and weights
% CFL condition (time-step restriction of RKDG schemes)
CFL = 0.7/(2*N+1);

%% Compute element mass matrix using Gaussian quadrature
Me = zeros(N+1,N+1); 
for j=1:nGP
    [phi,phi_xi]=BaseFunc(xiGP(j)); 
    for k = 1:N+1
        for l = 1:N+1
            Me(k,l) = Me(k,l) + wGP(j)*phi(k)*phi(l);
        end
    end
end
% Inverse of the element mass matrix
iMe = inv(Me);
%% Compute basis function for the grid used to plot
xisub = linspace(0,1,nSub);
phisub = zeros(N+1,nSub);
for i = 1:nSub
    [phi,phi_xi] = BaseFunc(xisub(i));
    phisub(:,i) = phi(:);
end
%% Store basis functions in Gaussian quadrature points 
% (frequently needed in the scheme)
for j = 1:nGP
    [phi,phi_xi] = BaseFunc(xiGP(j));
    phiGP(:,j) = phi(:);
    phi_xiGP(:,j) = phi_xi(:);
end
[phiL,phi_xi] = BaseFunc(0);
[phiR,phi_xi] = BaseFunc(1);

%% Compute initial condition and prepare fine output (subgrid 2N+1 points)
% Allocate degrees of freedom of the DG scheme
uhat = zeros(N+1,IMAX);
count = 0;
for i = 1:IMAX
    for j = 1:nGP
        xGP = x(i)-dx/2+xiGP(j)*dx; % Physical coordinate of the Gaussian point
        switch(ictype)
            case(0) % RP
                if (xGP<xc)
                    q0 = qL;
                else
                    q0 = qR;
                end
            case(1) % Function
                q0 = h(xGP);
        end
        [phi,phi_xi] = BaseFunc(xiGP(j)); 
        for k = 1:N+1
            uhat(k,i) = uhat(k,i) + wGP(j)*phi(k)*q0;
        end
    end
    uhat(:,i) = iMe*uhat(:,i);
    % Compute data to plot
    for j = 1:nSub
        count = count+1;
        xsub(count) = x(i)-dx/2+xisub(j)*dx; % Physical coordinate of subgrid point
        usub(count) = phisub(:,j)'*uhat(:,i);    % Initial solution
    end
end
%% Plot IC
% Plot only cell avearges
% plot(x,squeeze(uhat(1,:)),'o'); 
% Plot DG polynomials
figure
plot(xsub,usub,'b-') 

%% Initialize auxiliary variable to compute solution
uhatnew = uhat;

%% Start main time-loop
for n=1:NMAX
    %% Compute time step using CFL condition
    % Compute maximum signal speed
    amax = 0;
    for i=1:IMAX
        amax = max(amax,abs(a(uhat(1,i))));
    end
    % Compute time step
    dt = CFL*dx/amax;
    % Adjust the last time step to reach tend exactly
    if (time+dt>tend)
        dt = tend-time;
    end
    % Stop criterium
    if (time>=tend)
        break
    end
    
    %% Runge-Kutta scheme
    switch(rktype)
        case(1)
            % first order forward Euler scheme (unstable for high order DG schemes)
            k1 = LhDG(uhat);
            uhatnew = uhat + dt*k1; 
        case(2) 
            % second order RK scheme (unstable for high order DG schemes)
            k1 = LhDG(uhat);
            k2 = LhDG(uhat+dt/2*k1); 
            uhatnew = uhat + dt*k2; 
        case(3) 
            % 3rd order TVD Runge-Kutta scheme 
            uhat1   =     uhat                 + dt*LhDG(uhat);     
            uhat2   = 3/4*uhat + 1/4*uhat1 +   dt/4*LhDG(uhat1); 
            uhatnew = 1/3*uhat + 2/3*uhat2 + 2/3*dt*LhDG(uhat2);   
        case(4)
            % classical fourth order RK scheme
            k1 = LhDG(uhat); 
            k2 = LhDG(uhat+0.5*dt*k1); 
            k3 = LhDG(uhat+0.5*dt*k2); 
            k4 = LhDG(uhat+1.0*dt*k3);
            uhatnew = uhat + dt/6*(k1 + 2*k2 + 2*k3 + k4); 
        case(5)
            % use Butcher's RK5 scheme     
            k1 = LhDG(uhat); 
            k2 = LhDG(uhat+dt/4*k1);  
            k3 = LhDG(uhat+dt/8*k1+dt/8*k2);  
            k4 = LhDG(uhat-dt/2*k2+dt*k3); 
            k5 = LhDG(uhat+dt*3/16*k1+dt*9/16*k4); 
            k6 = LhDG(uhat-dt*3/7*k1+dt*2/7*k2+dt*12/7*k3-dt*12/7*k4+dt*8/7*k5); 
            uhatnew = uhat + dt/90*(7*k1+32*k3+12*k4+32*k5+7*k6);  
    end
    
    %% Overwrite solution
    uhat = uhatnew;
    %% Update time
    time = time +dt;
    
    %% Plot solution
    % Simply plot the cell averages
    % plot(x,squeeze(uhat(1,:)),'o') 
    % Plot DG polynomials
    count = 0;
    for i=1:IMAX
        for j = 1:nSub
            count = count+1;
            usub(count) = phisub(:,j)'*uhat(:,i);    % Initial solution
        end
    end
    plot(xsub,usub,'b-') 
    title(sprintf('Current time = %f',time))
    drawnow % Force Matlab to plot inmediately
end

%% Plot exact solution and compute errors
xe = linspace(xL,xR,10*IMAX);
switch(ictype)
    case(0)
        % calculate the exact solution on the fine mesh 
        for i=1:length(xe)
            qe(i) = ExactRiemannScalar(qL,qR,(xe(i)-xc)/time); 
        end
    case(1)
        % calculate the exact solution on the fine mesh 
        for i=1:length(xe)
            qe(i) = ExactCauchyScalar(xe(i),time); 
        end
        %% Calculate the L2 error norm 
        L2error = 0; 
        Linferror = 0; 
        for i=1:IMAX
            for j=1:nGP                 
                xGP = x(i)-dx/2+xiGP(j)*dx; % physical coordinate of the Gaussian quadrature point  
                qGP = phiGP(:,j)' * uhat(:,i); 
                qeGP = ExactCauchyScalar(xGP,time); 
                error = abs(qGP - qeGP); 
                L2error = L2error + dx*wGP(j)*error^2; 
                Linferror = max(Linferror, error); 
            end
        end
        L2error = sqrt(L2error); 
        disp(sprintf(' L2error   = %20.10e ', L2error)) 
        disp(sprintf(' Linferror = %20.10e ', Linferror)) 
end
hold on 
plot(xe,qe,'r-')
legend('RKDG scheme', 'Exact solution')


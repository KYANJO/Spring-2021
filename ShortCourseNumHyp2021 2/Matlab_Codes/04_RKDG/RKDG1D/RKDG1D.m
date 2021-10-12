%% Ruge-Kutta DG scheme for scalar conservation law
% q_t +f_x = 0
clear all
close all
clc

global N    % Polynomial degree of the basis functions
global IMAX nGP wGP xiGP Me iMe dx phiL phiR phiGP phi_xiGP

%% Define computational domain
xL = -1;
xR = 1;
N = 3; 
IMAX = 50;
time = 0;
tend = 0.5;
NMAX = 100000;

%% Select scheme
rktype = 3;

%% Define initial condition
% 0 == RP
% 1 == smooth function
ictype = 1;
qL = 1;
qR = 2;
xc = -0.5;

%% Define mesh, DOF ...
dx = (xR-xL)/IMAX;
x = linspace(xL+dx/2,xR-dx/2,IMAX);
nSub = 2*(N+1);  % Subgrgid for plotting
nGP = N+1;
nDOF = N+1;
[xiGP,wGP] = gauleg(nGP);
% CFL 
CFL = 0.7/(2*N+1);

% Compute elemnt mass matrix using Gaussian quadrature
Me = zeros(N+1,N+1);
for j=1:nGP
    [phi,phi_xi] = BaseFunc(xiGP(j));
    for k=1:N+1
        for l=1:N+1
            Me(k,l) = Me(k,l) +wGP(j)*phi(k)*phi(l);
        end
    end
end
% Inverse
iMe = inv(Me);

%% Compute basis functions on the subgrid
xisub = linspace(0,1,nSub);
phisub = zeros(N+1,nSub);
for i=1:nSub
    [phi,phi_xi] = BaseFunc(xisub(i));
    phisub(:,i) = phi(:);
end

%% Save basis functions in the Gaussian quadrature points
for j=1:nGP
    [phi,phi_xi] = BaseFunc(xiGP(j));
    phiGP(:,j) = phi(:);
    phi_xiGP(:,j) = phi_xi(:);
end
[phiL,phi_xi] = BaseFunc(0);
[phiR,phi_xi] = BaseFunc(1);

%% Compute the IC and prepare the final output on the subgrid (2N+1 points)
uhat = zeros(N+1,IMAX);
count = 0;
for i=1:IMAX
    for j=1:nGP
        xGP = x(i) -dx/2 + xiGP(j)*dx;
        switch(ictype)
            case(0) % RP
                if (xGP<xc)
                    q0=qL;
                else
                    q0 = qR;
                end
            case(1) % Smooth function
                q0 = h(xGP);
        end
        [phi,phi_xi] = BaseFunc(xiGP(j));
        for k = 1:N+1
            uhat(k,i) = uhat(k,i) + wGP(j)*phi(k)*q0;
        end
    end
    uhat(:,i) = iMe*uhat(:,i);
    % Compute data to plot
    for j=1:nSub
        count= count+1;
        xsub(count) = x(i) - dx/2 + xisub(j)*dx;
        usub(count) = phisub(:,j)'*uhat(:,i);
    end
end
% Plot IC
figure
plot(xsub,usub,'bo')

% Initialize auxiliary variable to compute the solution
uhatnew = uhat;

% Time loop
for n=1:NMAX
    %% Compute time step
    amax = 0;
    for i=1:IMAX
        amax = max(amax, abs(a(uhat(1,i))) );
    end
    dt = CFL*dx/amax;
    if (time +dt >tend)
        dt = tend-time;
    end
    if (time >= tend)
        break
    end
    
    %% Ruge-Kutta scheme
    switch(rktype)
        case(1) % (unstable for high order DG schemes)
        case(2) % (unstable for high order DG schemes)
        case(3)
            % 3d order RK scheme
            uhat1 = uhat +dt*LhDG(uhat);
            uhat2 = 3/4*uhat+1/4*uhat1+dt/4*LhDG(uhat1);
            uhatnew = 1/3*uhat+2/3*uhat2+2/3*dt*LhDG(uhat2);
        case(4)
        case(5)
    end
    % Overwrite solution
    uhat = uhatnew;
    % Update time
    time = time +dt;
    
    %% Plot solution
    count = 0;
    for i=1:IMAX
        for j=1:nSub
            count= count+1;
            usub(count) = phisub(:,j)'*uhat(:,i);
        end
    end
    plot(xsub,usub,'bo')
    title(sprintf('Time = %f',time))
    drawnow
end

%% Plot exact solution
xe = linspace(xL,xR,10*IMAX);
switch(ictype)
    case(0)
        
    case(1)
        for i=1:length(xe)
            qe(i) = ExactCauchyScalar(xe(i),time);
        end
        L2error = 0;
        Linf = 0;
        for i=1:IMAX
            for j=1:nGP
                xGP = x(i)-dx/2+xiGP(j)*dx;
                qeGP = ExactCauchyScalar(xGP,time);
                qGP = phiGP(:,j)'*uhat(:,i);
                error = abs(qGP-qeGP);
                L2error = L2error + dx*wGP(j)*error^2;
            end
        end
end
L2error= sqrt(L2error);
fprintf('L2error %e\n',L2error)
hold on
plot(xe,qe,'b-')

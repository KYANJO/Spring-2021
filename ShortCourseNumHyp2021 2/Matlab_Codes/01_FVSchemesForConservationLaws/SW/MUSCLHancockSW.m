% FV solver, based on MUSCL-Hancock
clear all
close all
clc

global g    
g = 9.81;   % Gravity constant

%% Define problem
% hL = 1;   % Left water depth
% hR = 0.2; % Right water depth
% uL = 1;   % Left velocity
% uR = -1;  % Right velocity

hL = 2;     % Water depth a sinistra
hR = 0.001; % Water depth a destra
uL = 0;     % Left velocity
uR = 0;     % Right velocity

psiL = 1;
psiR = 0.5;
QL = [hL; hL*uL; hL*psiL]; % Left conservative vector
QR = [hR; hR*uR; hR*psiR]; % Right conservative vector

time = 0;       % Initial time
tend = 0.1;     % End time
xL = -1; xR=1;  % Computational domain

IMAX = 200;     % Number of control volues
dx = (xR-xL)/IMAX; 
x = linspace(xL+dx/2,xR-dx/2,IMAX);
CFL = 0.5;
NMAX = 10000;

%% Define initial condition
Q = zeros(3,IMAX);
for i=1:IMAX
    if (x(i)<=0) 
        Q(:,i) = QL;
    else
        Q(:,i) = QR;
    end
end
subplot(3,1,1)
plot(x,Q(1,:))             % h
subplot(3,1,2)
plot(x,Q(2,:)./Q(1,:))    % u
subplot(3,1,3)
plot(x,Q(3,:)./Q(1,:))    % psi

%% Compute approximate solution
for n = 1:NMAX % loop in tempo
    % Compute the maximum of the abs eigenvalues
    amax = 0;
    for i=1:IMAX
        L = Lambda(Q(:,i));
        amax = max(amax,max(abs(L)));
    end
    % Compute time step
    dt = CFL*dx/amax;    
    % Controlliamo di arrivare al tempo finale essatamente
    if (time+dt>tend)
        dt = tend - time;
    end
    % Stop criterium
    if (time>=tend)
        break
    end
    
    %Piecewise linear reconstuction in space
    slope = zeros(3,IMAX); WR = zeros(3,IMAX); WL = zeros(3,IMAX);
    for i=1:IMAX
        if (i==1)
            slope(:,i) = zeros(3,1);
        elseif (i==IMAX)
            slope(:,i) = zeros(3,1);
        else % Minmod
            slope(:,i)  = minmod(Q(:,i)-Q(:,i-1),Q(:,i+1)-Q(:,i));
        end  
        % Compute extrapolated values at the cell boundary
        WR(:,i) = Q(:,i) + 0.5.*slope(:,i); % valore a destra
        WL(:,i) = Q(:,i) - 0.5.*slope(:,i); % valore a sinistra
        % Compute time derivative (Cauchy-Kovalevskaya): Q_t=-f_x
        Qt = - (f(WR(:,i))-f(WL(:,i)))/dx;
        % Update extrapolated values with the half in time evolution
        WR(:,i) = WR(:,i) + 0.5*dt*Qt;
        WL(:,i) = WL(:,i) + 0.5*dt*Qt;
    end
    
    % Loop on control volumes
    for i=1:IMAX
        if (i==1)
             % Dirichlet
             %Qnew(:,i) = QL;
             % Reflective BC
             QGod = ExactRiemannSW(WR(:,i),WL(:,i+1),0);
             fp = f(QGod);
             QBC = WL(:,i); % copy the inner state
             QBC(2) = - QBC(2); % invertiamo la velocita normale
             QBC(3) = - QBC(3);
             QGod = ExactRiemannSW(QBC,WL(:,i),0);
             fm = f(QGod);
             Qnew(:,i) = Q(:,i) - dt/dx*(fp-fm);
        elseif (i==IMAX)
             % Dirichlet
             %Qnew(:,i) = QR;
             % reflective
             QBC = WR(:,i); % copy the inner state
             QBC(2) = - QBC(2); % invertiamo la velocita normale
             QBC(3) = - QBC(3);
             QGod = ExactRiemannSW(WR(:,i),QBC,0);
             fp = f(QGod);
             QGod = ExactRiemannSW(WR(:,i-1),WL(:,i),0);
             fm = f(QGod);
             Qnew(:,i) = Q(:,i) - dt/dx*(fp-fm);
        else
            % Godunov flux
            QGod = ExactRiemannSW(WR(:,i),WL(:,i+1),0);
            fp = f(QGod);
            QGod = ExactRiemannSW(WR(:,i-1),WL(:,i),0);
            fm = f(QGod);
            % FV method
            Qnew(:,i) = Q(:,i) - dt/dx*(fp-fm);
        end
    end  
    % Update time
    time = time + dt;
    % Overrite solution
    Q = Qnew;
    
    % Plot
    subplot(3,1,1)
    plot(x,Q(1,:),'bo')             % h
    title(sprintf('Time = % f',time))
    subplot(3,1,2)
    plot(x,Q(2,:)./Q(1,:),'bo')    % u
    subplot(3,1,3)
    plot(x,Q(3,:)./Q(1,:),'bo')    % psi
    drawnow
end

%% Plot exact solution
% Compute solution at tend with the exact Riemann solver
for i= 1:IMAX
    xi = x(i)/time;
    Qe(:,i) = ExactRiemannSW(QL,QR,xi);
end
% Plot the exact solution at the final time
hold all
subplot(3,1,1)
hold on
plot(x,Qe(1,:),'k-')             % h
subplot(3,1,2)
hold on
plot(x,Qe(2,:)./Qe(1,:),'k-')    % u
subplot(3,1,3)
hold on
plot(x,Qe(3,:)./Qe(1,:),'k-')    % psi
drawnow
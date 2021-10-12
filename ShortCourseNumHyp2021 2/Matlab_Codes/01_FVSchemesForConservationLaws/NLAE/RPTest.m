close all
clear all
clc 
%% Program  for the visualization of the exact solution of the 
% Riemann problem for the scalar equation
% q_t + f_x = 0

%Definition of the physical problem
xL = -1; 
xR = 1;
time =  0; 
tend = 0.4;
dt = 1.e-3;
% Number of spatial nodes
IMAX = 100;
% Discretization in space
 x = linspace(xL,xR,IMAX);
% Number of time steps
NMAX = tend/dt;
% Update of the time
time = time + dt;

% Initial condition
qL = 2; qR = 1; % Initial states
%qL = -0.5; qR = 2; % Initial states

%% Time loop
for n = 1:NMAX
    % Loop in nodes
    for i = 1:IMAX
        % Compute similarity parameter
        xi = x(i)/time;
        % Compute exact solution
        q(i) = ExactRiemannSolverNL(qL,qR,xi);
    end
    % Plot
    plot(x,q,'r-')
    title(sprintf('Time = %f',time))
    drawnow
    % Time update
    time = time + dt;
end
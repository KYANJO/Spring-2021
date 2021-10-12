clear all
close all
%figure
clc
global g    % gravity constant
g = 9.81;   % 
% dambreak problem 
hL = 1    % left depth
hR = 0.1      % right depth
uL = 0   % left velocity
uR = 0    % right veocity
psiL = 1;   % left concentration
psiR = 0;   % right concentration
hs = Newton(hL,hR,uL,uR) 


% plot the exact solution of the Riemann problem 
QL = [hL;hL*uL;hL*psiL];  % left conservative variables 
QR = [hR;hR*uR;hR*psiR];  % right conservative variables 
xL = -1; 
xR = +1; 
time = 0.1; 
IMAX = 1000; 
x = linspace(xL,xR,IMAX); 
for i=1:IMAX
    xi = x(i)/time; 
    Qe(:,i) = ExactRiemannSW(QL,QR,xi);    
end
subplot(3,1,1) 
plot(x,Qe(1,:))          % h 
subplot(3,1,2) 
plot(x,Qe(2,:)./Qe(1,:)) % u 
subplot(3,1,3) 
plot(x,Qe(3,:)./Qe(1,:)) % psi 
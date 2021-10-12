%Sylvester equations
clc
close all

SA = [4 7 -6 10 9; 4 -6 4 9 5; -2 4 6 10 3; -4 6 3 -3 7;...
       -1 8 0 6 2];
va = [  4, 1,3,9,10]; VA = diag(va);

SB = [6 6 -1 5;8 7 -6 -6; -3 3 -5 10; -6 -6 -9 -7];

vb = [-7, -4, -3, -5]; VB = diag(vb);

C = [-9 10 6 -7; -8 -2 -5 3; -7 0 -6 5;-8 9 0 8;-4 -9 -4 5];


%format long
[m,n] = size(C);

%identity matrix
e = ones(m,1); I = diag(e);

%inverse of SB and SA
invSB = inv(SB); invSA = inv(SA);

xhat = zeros(m,n); x = zeros(m,n);

for i = 1:n
   A = VA - VB(i,i)*I; 
   Chat = SA*C*invSB(:,i);
   xhat(:,i) = A\Chat;
end

J = 1:n;
fprintf("The solution in matrix form is:\n\n")
x(:,J) = invSA*xhat*SB(:,J)


   

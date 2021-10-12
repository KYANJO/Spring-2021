%Sylvester equations
clear all
close all
clc

SA = [4 7 -6 10 9; 4 -6 4 9 5; -2 4 6 10 3; -4 6 3 -3 7;-1 8 0 6 2];

inv_SA = inv(SA);

va = [4,1,3,9,10]; VA = diag(va);

SB = [6 6 -1 5;8 7 -6 -6; -3 3 -5 10; -6 -6 -9 -7];
inv_SB = inv(SB);

vb = [-7, -4, -3, -5]; VB = diag(vb);

C = [-9 10 6 -7; -8 -2 -5 3; -7 0 -6 5;-8 9 0 8;-4 -9 -4 5];

[m,n] = size(C);
I = diag(ones(m,1));
Xhat = zeros(m,n);

for j = 1:n

 M = VA - VB(j,j)*I;
 b = SA*C*inv_SB(:,j);
 Xhat(:,j) = M\b;

end


X = zeros(m,n);

for j = 1:n
      
      X(:,j) = inv_SA*Xhat*SB(:,j)
      
end

A = inv_SA*VA*SA;
B = inv_SB*VB*SB;

X1 = sylvester(A,-B,C);
 
 er = X-X1;
 
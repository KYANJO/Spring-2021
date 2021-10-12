clear all
close all
clc
% Dimension of the matrix
m = 50;   
n = 12;

% Vendermond matrix
A = zeros(m,n);

J = 1:m;
t = (J-1)/(m-1);

for j = 1:n
    
    A(:,j) = t.^(j-1);
end

f = @(t) cos(4.*t);

b = f(t)';

% use of the new method for solving the least square problem
I = eye(m);
Z = zeros(n);
B = [b; zeros(n,1)];
AA = [I A; A' Z];
% compute solution xr = [r ; x]
%xr = AA\B;
[Q,R] = qr(AA);
y_q = Q'*B;
xr = R\y_q;
% extract the solution for the least square problem from xr
x = xr(m+1:end);
r = xr(1:m);
% (a) Use of the normal equations
x_normal = (A'*A)\(A'*b);
% (e) QR builtin function
[Q1,R1] = qr(A);
y = Q1'*b;
x_qr = R1\y;

% Table of coefficients for each methods

Tab = table(round(x,16),round(x_normal,16),round(x_qr,16),'VariableNames',...
    {'New method','Normal(a)','QR Matlab(e)'});

disp(Tab)




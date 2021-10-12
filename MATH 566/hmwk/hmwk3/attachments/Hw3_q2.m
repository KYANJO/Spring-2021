clear all
close all
clc
% Dimension on the matrix
m = 50;   
n = 12;
gamma = 1;
% Vendermond matrix
A = zeros(m,n);
I = eye(n);
J = 1:m;
t = (J-1)/(m-1);

for j = 1:n
    
    A(:,j) = t.^(j-1);
end

gw = @(t, tj, gamma) exp(-(abs(t-tj)/gamma).^2);
f = @(t) cos(4.*t);

b = f(t)';

% QR builtin function
[Q,R] = qr(A);
y = Q'*b;
% non-weighted solution
x = R\y;

tf = 1/23;

W = diag(gw(tf,t, gamma));

[Q1,R1] = qr(W*A);
y1 = Q1'*(W*b);
% weighted solution
x1 = R1\y1;
fprintf('Coefficients of the weighted least squares solution x:\n\n');
disp(x1')
% p11(1/23)
K = 1:n;
tp = tf.^(K-1);
p1 = sum(x1'.*tp);
%fprintf('Value of the polynomial with weight at t= 1/23\n\n');
%disp(round(p1,15))

p = sum(x'.*tp);
%fprintf('Value of the polynomial with non-weight at t= 1/23\n\n');
%disp(round(p,15))

fprintf('Value of the polynomials at t= 1/23\n\n');
Tab = table(1/23, round(p,15),round(p1,15),'VariableNames',{'t',...
    'Non-weighted','weighted'});

disp(Tab)

fprintf('The value of the non-weighted polynomial at t = 1/23 is less\n');
fprintf('compared to that of weighted polyomial.\n\n');
fprintf('Weighted least squares method provides a better approximation \n');
fprintf('to f(t) at t= 1/23.\n');
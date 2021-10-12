
% Compute the roots of the sixth degree Legendre polynomial
% p(z) = (1/16)*(231*z.^6 - 315*z.^4 + 105*z.^2 -5);
clear all
close all
clc

%functions
q = @(z) (1/16)*(231*z.^6 - 315*z.^4 + 105*z.^2 -5);
%monic polynomial
p = @(z) (16/231)*q(z);

%companion matrix
syms z

% coefficients of the p(z)
coef = coeffs(p(z),'all');
% degree of the polynomial
n = polynomialDegree(p(z));
% rearangement of the coefficients to compute the companion matrix
coef = fliplr(coef);
coef = coef(1:n);
coef = double(coef);
% companion matrix
I = eye(n-1,n-1); 
A = [zeros(n-1,1) I];
A = [A; -coef];
% roots of the polynomial p(z)
roots = eig(A);
% verify that those roots are the actual roots of the polynomial
ql = q(roots)

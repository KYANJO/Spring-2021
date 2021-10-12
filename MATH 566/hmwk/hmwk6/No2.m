%Companion matrix
% Compute the roots of the sixth degree Legendre polynomial
clc
close all

fprintf("No.2\n\n");

%functions
q = @(z) (1/16)*(231*z.^6 - 315*z.^4 + 105*z.^2 -5);
%monic polynomial
p = @(z) (16/231)*q(z);

%companion matrix
syms z
A = companion(p(z));

%eigen values
format long
roots = eig(A);
fprintf("Roots in ascending order:\n");
roots = sort(roots);
disp(roots);
fprintf("Check if they are actual roots\n");
q(roots)
fprintf("Hence the roots are actual\n");

function A = companion(p)
    %coeffcients of poly
    C = coeffs(p,'all');
    cof = fliplr(C);
    %degree of poly
    n = polynomialDegree(p);
    cof = (cof(1:n));
    cof = double(cof);
    I = eye(n-1,n-1); 
    A = [zeros(n-1,1) I];
    A = [A; -cof];

end

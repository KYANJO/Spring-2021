function x = jacobi(A, b, x0, steps, omega)
% Jacobi smoother
%
% INPUT:   A      matrix
%          b      right hand side
%          x0     initial guess
%          steps  number of iteration steps
%          omega  relaxation parameter
%
% OUTPUT:  x      approximate solution
%
% VERSION 1.0
% DATE 25.3.2004
% EMAIL bernd@flemisch.net

  d = diag(A);
  dinv = 1./d;
  x = x0;
  for k = 1:steps
    x = x + omega*dinv.*(b - A*x);
  end
  
  
  
  
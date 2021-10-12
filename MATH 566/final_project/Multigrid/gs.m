function x = gs(A, b, x0, steps, omega)
% SOR smoother
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
  
  n = size(A, 1);
  if omega == 1 
    N = tril(A);
    M = -triu(A,1);
  else
    D = spdiags(spdiags(A,0),0,n,n);
    L = tril(A,-1);
    N = D + omega*L;
    M = (1 - omega)*D - omega*triu(A,1);
    b = omega*b;
  end
  x = x0;
  for k = 1:steps
    x = N\(b + M*x);
  end
  
  
  
  
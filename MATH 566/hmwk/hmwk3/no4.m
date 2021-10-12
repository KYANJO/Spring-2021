%Stability of linear systems.

fprintf('4a).\n');
%(a).Does this imply that xˆ is close to the exact solution x?
fprintf('No, since the residual values are not near zero, then xhat is not close to the exact solution. \n\n');

fprintf('4b).\n');
%Matrix A
A = 1./hankel(2:6,6:10);

%vector b
b = [0.882 0.744 0.618 0.521 0.447]';

%Accurate solution to the system
x = A\b

fprintf('\n'); 
fprintf('4c).\n');
%Obtain a condition number for A using this same software again
%condition number 
C = cond(A)
fprintf('Since the Condition number of A is large then the system is ill-conditioned, therefore a small perturbation \n to the RHS can lead to large change ib the system. \n\n')
 
fprintf('Consider a small perturbation on, db. \n');

db = [0.000002 0.000004 0.000008 0.00001 0.00007]'

R1=C*norm(db,2)/norm(b,2)

%ddx due to perturbation on the RHS
ddx = x + A\db

%Relative error
RE = norm((ddx-x),2)/norm(x,2)

fprintf('Since the relative Error,RE<=R1, and large then indeed this confirms that a very small residual\n after the system being perturbed on the RHS, is big enough to allow for the solution to be as far away. \n')

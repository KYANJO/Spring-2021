function [u_lmax, k, flag] = mg_alg(f_lmax, lmax, matrices, u0_lmax, parents, ...
				       isdir, nu1, nu2, mu, maxit, tol, smooth, omega)
% Multigrid algorithm. 
% Can be used as standalone solver or as preconditioner for the
% Matlab pcg routine (see the pcg help for details).
%
% INPUT:   f_lmax    right hand side on the finest level
%          lmax      maximum level
%          matrices  cell array of matrices, matrices{l+1} contains
%                    the system matrix on level l, l = 0,...,lmax
%          u0_lmax   initial guess
%          parents   cell array containing the father-son relations
%          isdir     dirichlet node flags
%          nu1       number of pre-smoothing steps
%          nu2       number of post-smoothing steps
%          mu        number of recursive multigrid step calls, 
%                    1 for V-cycle, 2 for W-cycle
%          maxit     maximum number of iterations
%          tol       tolerance for stopping criterion
%          smooth    type of smoother, 
%                    'Jacobi' for damped Jacobi or 'SOR' for SOR
%          omega     relaxation parameter for the smoother
%
% OUTPUT:  u_lmax    solution on the finest level
%          k         number of performed iterations
%          flag      0 if the algorithm converged, 1 otherwise
%
% VERSION 1.0
% DATE 25.3.2004
% EMAIL bernd@flemisch.net

  flag = 0;
  for k = 1:maxit
    % perform a multigrid step on level lmax
    u_lmax = multigrid_step(f_lmax, lmax, matrices, u0_lmax, parents, ...
			 isdir, nu1, nu2, mu, smooth, omega);
    % stopping criterion (relative residual):
    if (norm(matrices{lmax+1}*u_lmax - f_lmax)/norm(f_lmax) < tol) 
      return;
    end
    u0_lmax = u_lmax;
  end
  flag = 1;
  return;
  
function u_l = multigrid_step(f_l, l, matrices, u0_l, parents, isdir, ...
			      nu1, nu2, mu, smooth, omega)
  A_l = matrices{l+1}; % get system matrix of the current level
  if (strcmp(smooth, 'Jacobi'))
    u_l = jacobi(A_l, f_l, u0_l, nu1, omega); % smooth nu1-times
  else 
    u_l = gs(A_l, f_l, u0_l, nu1, omega); % smooth nu1-times
  end
  d_l = f_l - A_l*u_l; % calculate defect
  d_lm1 = restrict(d_l, parents{l+1}, isdir); % restrict defect
  if l == 1     % solve direct on the coarsest level:
    w_lm1 = matrices{l}\d_lm1;
  else          % perform mu multigrid steps
    u0_lm1 = zeros(length(d_lm1), 1);
    for k = 1:mu
      w_lm1 = multigrid_step(d_lm1, l-1, matrices, u0_lm1, parents, ...
			     isdir, nu1, nu2, mu, smooth, omega);
      u0_lm1 = w_lm1;
    end
  end
  w_l = interpolate(w_lm1, parents{l+1}); % interpolate correction
  u_l = u_l + w_l; % add correction
  if (strcmp(smooth, 'Jacobi'))
    u_l = jacobi(A_l, f_l, u_l, nu2, omega); % smooth nu2-times
  else 
    u_l = gs(A_l, f_l, u_l, nu2, omega); % smooth nu2-times
  end
    
  
  
  
  
  
function [p, e, t, u_lmax] = mg_main(steps, Hmax, smooth, omega, nu1, nu2, mu, tol, ...
				     maxit, exact, direct, prec, geom, b, a, c, f, u)
% Main routine for multigrid solver. 
% Called by mg_menu.m and mg_test.m
% 
% INPUT:   steps   number of refinement steps
%          Hmax    maximum initial element diameter
%          smooth  type of smoother, 
%                  'Jacobi' for damped Jacobi or 'SOR' for SOR
%          omega   relaxation parameter for the smoother
%          nu1     number of pre-smoothing steps
%          nu2     number of post-smoothing steps
%          mu      number of recursive multigrid step calls, 
%                  1 for V-cycle, 2 for W-cycle
%          tol     tolerance for stopping criterion
%          maxit   maximum number of multigrid iterations
%          exact   flag whether exact solution is available
%          direct  flag whether a direct solve is performed
%          prec    flag whether multigrid is used as preconditioner
%          geom    Matlab geometry description
%          b       Matlab boundary condition description
%          a       problem parameters for the PDE toolbox to solve
%          c       div(a*grad(u)) + c*u = f
%          f       right hand side
%          u       exact solution
%
% OUTPUT:  p,e,t   mesh data on the finest level
%          u_lmax  FE solution vector on the finest level
%
% VERSION 1.0
% DATE 25.3.2004
% EMAIL bernd@flemisch.net

matrices = cell(steps+1,1); % matrices{l+1} contains the system
                            % matrix on level l, l = 0,...,lmax
parents = cell(steps+1,1); % parents{l+1} contains the father-son information
			   % on level l, needed for the transfer operators
time_dir = 0; error_dir = 0; error = 0;
[p, e, t] = initmesh(geom, 'Hmax', Hmax); % generate initial mesh
dof = size(p, 2); % number of degrees of freedom
[matrices{1}, f_lmax] = assempde(b, p, e, t, a, c, f); % assemble system on level 0
if (prec) fprintf('MG-PCG: '); else fprintf('Multigrid: '); end
if (mu == 1) fprintf('V-cycle, '); else fprintf('W-cycle, '); end
if strcmp(smooth, 'Jacobi') fprintf('Jacobi-'); 
else fprintf('SOR-'); end
fprintf('Smoother (%.2f), pre = %d, post = %d, tol = %.1e\n', ...
	omega, nu1, nu2, tol);
fprintf('lmax \t dof\t iter \t error MG\t error direct\t time MG\t time direct\n');

for lmax = 1:steps
  isdir = dirichlet_flags(p, e, b); % get dirichlet flags
  parents{lmax+1} = get_edges(p, t);  % get the father-son information
  [p, e, t] = refinemesh(geom, p, e, t); % refine mesh uniformly
  dof = size(p, 2); % number of degrees of freedom
  [matrices{lmax+1}, f_lmax] = assempde(b, p, e, t, a, c, f); 
                                      % assemble system on level lmax
  u0_lmax = zeros(dof,1); % initial guess
  if (prec) % use multigrid as CG preconditioner, see the pcg help
            % for details
    tic;
    [u_lmax, flag, rr, iter] = pcg(matrices{lmax+1}, f_lmax ,tol, maxit, @mg_alg, ...
				   [],[], lmax, matrices, u0_lmax, parents, ...
				   isdir, nu1, nu2, mu, 1, eps, smooth, omega);
    time = toc;
  else      % use multigrid as solver
    tic;
    [u_lmax, iter, flag] = mg_alg(f_lmax, lmax, matrices, u0_lmax, parents, isdir, ...
				  nu1, nu2, mu, maxit, tol, smooth, omega);
    time = toc;
  end
  if (flag)
    fprintf('Multigrid did not converge within %d iterations.\n', maxit);
  end
  if (direct)
    tic;
    u_lmax_dir = matrices{lmax+1}\f_lmax; % direct solve
    time_dir = toc;
  end
  if (exact)
    x = p(1,:);
    y = p(2,:);
    exactsol = eval(u);
    %exactsol = feval(u, p(1,:), p(2,:)); % evaluate exact solution
    error = norm(u_lmax - exactsol'); % calculate the error in the
                                      % euclidean vector norm
    if (direct)
      error_dir = norm(u_lmax_dir - exactsol');
    end
  end
  fprintf('%d \t %d\t %d \t %8.2e \t %8.2e \t %8.2e \t %8.2e\n', ...
	  lmax, dof, iter, error, error_dir, time, time_dir);
end








% Numerical approximation to Poisson's equation over the square [a,b]x[a,b] with
% Dirichlet boundary conditions.  Uses a uniform mesh with (n+2)x(n+2) total
% points (i.e, n interior grid points). 
%
% Solves with classical multigrid (v-cycle and damped-Jacobi smoother).
%
% Input:
%     ffun : the RHS of poisson equation (i.e. the Laplacian of u).
%     gfun : the boundary function representing the Dirichlet B.C.
%      a,b : the interval defining the square
%        m : m+2 is the number of points in either direction of the mesh.
% smoother : Selection of the smoother.
% Ouput:
%        u : the numerical solution of Poisson equation at the mesh points.
%      x,y : the uniform mesh.
%
function [u,x,y] = fd2poissonmg(ffun,gfun,a,b,m,smoother)

k = round(log2(m+1));

if abs(log2(m+1)-k) > 100*eps
    error('fd2poissonMg:power2','m+1 must be a power of 2');
end

h = (b-a)/2^k;   % Mesh spacing

[x,y] = meshgrid(a:h:b);   % Uniform mesh, including boundary points.

idx = 2:m+1;
idy = 2:m+1;

% Compute boundary terms, south, north, east, west
ubs = feval(gfun,x(1,1:m+2),y(1,1:m+2));     % Include corners
ubn = feval(gfun,x(m+2,1:m+2),y(m+2,1:m+2)); % Include corners
ube = feval(gfun,x(idy,m+2),y(idy,m+2));     % No corners
ubw = feval(gfun,x(idy,1),y(idy,1));         % No corners

% Evaluate the RHS of Poisson's equation at the interior points.
f = feval(ffun,x,y);

% Inital guess
u = zeros(m+2,m+2);

% Add boundary terms to u
u(idx,1) = ubw; u(idx,m+2) = ube;
u(1,:) = ubs; u(m+2,:) = ubn;

tol = 1e-8;
r = inf;
nrmf = norm(f(:));
while norm(r(:)) > tol*nrmf
    [u,r] = vcycle(u,f,h,k,smoother);
end

end

function [u,r] = vcycle(u,f,h,k,smoother)

numSmooths = 3;
n = 2^k; n2 = n/2;
i = 2:n; i2 = 2:n2;

if k == 1
    % Only one interior point so solve the system:
    u(2,2) = 0.25*(u(1,2)+u(3,2)+u(2,1)+u(2,3)-h^2*f(2,2));
    return;
else
    % Smooth the solution.
    %dampedJacobi
    %u = dampedJacobi(u,f,h,k,numSmooths);
    %Red-Black Gaus
    %u = redblack(u,f,h,k,numSmooths);

    if smoother == 'dampedJacobi'
        u = dampedJacobi(u,f,h,k,numSmooths);
    elseif  smoother == 'redblack'
        u = redblack(u,f,h,k,numSmooths);
    end
    
    % Calculate residual
    r = zeros(n+1);
    r(i,i) = f(i,i) - (-4*u(i,i)+u(i-1,i)+u(i+1,i)+u(i,i-1)+u(i,i+1))/h^2;
    
    % Restrict the residual to the coarser grid using full-weighting
    jj = 3:2:n-1;
    f2 = zeros(n2+1);
    f2(i2,i2) = (r(jj-1,jj)+r(jj+1,jj)+r(jj,jj-1)+r(jj,jj+1)+4*r(jj,jj))/8;
    
    % Repeat call to vcylce with the residual on coarser grid and with a
    % zero initial guess.
    u2 = vcycle(zeros(n2+1),f2,2*h,k-1);
    
    % Interpolate (prolongate) the correction to the finer grid
    ut = zeros(n+1);
    ut(1:2:n+1,1:2:n+1) = u2;
    ut(2:2:n,1:2:n+1) = 0.5*(u2(1:n2,:)+u2(2:n2+1,:));
    ut(:,2:2:n+1) = 0.5*(ut(:,1:2:n-1)+ut(:,3:2:n+1));
    
    % Update the solution and smooth again:
    u = dampedJacobi(u+ut,f,h,k,numSmooths);
    % Calculate residual
    r(i,i) = f(i,i) - (-4*u(i,i)+u(i-1,i)+u(i+1,i)+u(i,i-1)+u(i,i+1))/h^2;
end
end

%Damped Jacobi smoother
function u = dampedJacobi(u,f,h,k,numSmooths)
omega = 4/5;
n = 2^k;
i = 2:n;
f = h^2*f;
% Smooth u using underrelaxed Jacobi
for j=1:numSmooths
    u(i,i) = (0.25*omega)*(u(i-1,i)+u(i+1,i)+u(i,i-1)+u(i,i+1)-f(i,i)) + ...
                     (1-omega)*u(i,i);
end

end

%Red black Gauss-Seidel step
function u = redblack(u,f,h,k,numSmooths)
n = 2^k;
f = h^2*f;
for k = 1:numSmooths
    %Red points
    for j = 2:n+1
        for i = 2:n+1
          if (mod(i+j,2) == 0)
             u(i,j) = 1\4*(u(i-1,j)+u(i+1,j)+...
                 u(i,j-1)+u(i,j+1)+f(i,j));
          end
       end
    end

    %Black points
    for j = 2:n+1
       for i = 2:n+1
          if (mod(i+j,2) == 1)
             u(i,j) = 1\4*(u(i-1,j)+u(i+1,j)+...
                 u(i,j-1)+u(i,j+1)+f(i,j));
          end
       end
    end
end
end
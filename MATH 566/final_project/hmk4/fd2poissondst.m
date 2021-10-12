% Numerical approximation to Poisson's equation over the square [a,b]x[a,b] with
% Dirichlet boundary conditions.  Uses a uniform mesh with (n+2)x(n+2) total
% points (i.e, n interior grid points). 
%
% Solves with the DST.
%
% Input:
%     pfun : the RHS of poisson equation (i.e. the Laplacian of u).
%     bfun : the boundary function representing the Dirichlet B.C.
%      a,b : the interval defining the square
%        m : m+2 is the number of points in either direction of the mesh.
% Ouput:
%        u : the numerical solution of Poisson equation at the mesh points.
%      x,y : the uniform mesh.
%
function [u,x,y] = fd2poissondst(pfun,bfun,a,b,m)

h = (b-a)/(m+1);   % Mesh spacing

[x,y] = meshgrid(a:h:b);   % Uniform mesh, including boundary points.

idx = 2:m+1;
idy = 2:m+1;

% Compute boundary terms, south, north, east, west
ubs = feval(bfun,x(1,1:m+2),y(1,1:m+2));     % Include corners
ubn = feval(bfun,x(m+2,1:m+2),y(m+2,1:m+2)); % Include corners
ube = feval(bfun,x(idy,m+2),y(idy,m+2));     % No corners
ubw = feval(bfun,x(idy,1),y(idy,1));         % No corners

% Evaluate the RHS of Poisson's equation at the interior points.
f = feval(pfun,x(idy,idx),y(idy,idx));

% Adjust f for boundary terms
f(:,1) = f(:,1) - ubw/h^2;             % West
f(:,m) = f(:,m) - ube/h^2;             % East
f(1,1:m) = f(1,1:m) - ubs(idx)/h^2;    % South
f(m,1:m) = f(m,1:m) - ubn(idx)/h^2;    % North

% Computation of fhat=(S*f)*S^(-1), where S is the discrete sine transform.
fhat = idst(dst(f,1),2);

% Denominator for the computation of uhat:
denom = bsxfun(@plus,cos(pi*(idx-1)./(m+1)).',cos(pi*(idx-1)./(m+1)))-2;

uhat = h^2/2*(fhat./denom);

% Computation of u = (S^(-1)*uhat)*S
u = dst(idst(uhat,1),2);

% Append on to u the boundary values from the Dirichlet condition.
u = [ubs;[ubw,u,ube];ubn];

end




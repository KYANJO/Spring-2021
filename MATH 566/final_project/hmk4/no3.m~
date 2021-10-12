% Numerical approximation to poisson's equation over the square [a,b] x
% [a,b] with zero Neumann boundary conditions. Uses a unifoorm mesh with
% (n+2) x (n+2) total points.

% Solves with the DCT

% Input
%   pfun : the RHS of poisson equation (i.e. the Laplacian of u). (f(x,y))
%   bfun : the boundary function representing the Neumann B.C. 
%    a,b : the interval defining the square
%    m : m+2 is the number of points in either direction of the mesh.
%   

%Output
%   u : the numerical solution of poisson equation at the mesh points.
% x,y : the uniform mesh

function [u,x,y]=fd2poissondct(pfun,bfun,a,b,m)
h=(b-a)/(m+1); %mesh spacing
[x,y] = meshgrid(a:h:b); %uniform mesh, including boundary points.

idx = 2:m+1;
idy = 2:m+1;

%Compute boundary terms, south, north, east, west
ubs = 0 % feval(bfun, x(1,1:m+2),y(1,1:m+2));   %Include corners
ubn = 0 %feval(bfun, x(m+2,1:m+2),y(m+2,1:m+2)); %Include corners
ube = 0 %feval(bfun, x(idy,m+2),y(idy,m+2));       %no corners
ubw = 0 %feval(bfun, x(idy,1),y(idy,1));       %no corners

%Evaluate the RHS of Poisson's equation at the interior points.
f = feval(pfun, x(idy,idx),y(idy,idx));

%Adjust f for boundary terms

%Computation of fhat
fhat=idct(dct(f,1),2);

%Denominator for the computation of uhat:
demon = bsxfun(@plus,cos(pi*(idx-1)./(m+1).',cos(pi*(idx-1)./(m+1)))-2;

uhat = h^2/2*(fhat./demon);

%computation of u
u=dct(idct(uhat,1),2);

%Append on to u the boundary values from the Neumann condition.
u=[ubs;[ubw,u,ube];ubn];

end







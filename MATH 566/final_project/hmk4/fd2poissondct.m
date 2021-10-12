% Numerical approximation to poisson's equation over the square [a,b] x
% [a,b] with zero Neumann boundary conditions. Uses a unifoorm mesh with
% (n+2) x (n+2) total points.

% Solves with the DCT

% Input
%   pfun : the RHS of poisson equation (i.e. the Laplacian of u). (f(x,y))
%    a,b : the interval defining the square
%    m : m+2 is the number of points in either direction of the mesh.
%   

%Output
%   u : the numerical solution of poisson equation at the mesh points.
% x,y : the uniform mesh

function [u,x,y] = fd2poissondct(p,a,b,m)

h=1/(m+1);

% idx and idy need to include all the grid points:
idx = 1:m+2;
idy = 1:m+2;

[x,y] = meshgrid(a:h:b); %uniform mesh, including boundary points.

% Evaluate the RHS of Poisson's equation at the interior points.
fr = feval(p,x(idy,idx),y(idy,idx));

% Computation of fhat=(S*f)*S^(-1), where S is the DCT
fhat=idct(dct(fr,1),2);

% Denominator for the computation of uhat:
denom = [bsxfun(@plus,cos(pi*(idx-1)./(m+1)).',cos(pi*(idx-1)./(m+1)))-2];

uhat = h^2/2*(fhat./denom);

%Dealing with the zero eigenvalue.
uhat(1)=0;

% Computation of u = (S^(-1)*uhat)*S
u = dct(idct(uhat,1),2);

end
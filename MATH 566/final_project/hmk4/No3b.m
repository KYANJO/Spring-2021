% Solves the Poisson equation with f(x,y) = -8*(pi^2)*cos(2*pi*x)cos(2*pi*y)

m=(2^6)-1;
a=0;b=1;

%fuctionf(x,y)
p=@(x,y) -8*(pi^2)*(cos(2*pi*x)).*(cos(2*pi*y));

%Approximated
[u1,x,y]=fd2poissondct(p,a,b,m)

%Exact function
uex=@(x,y) (cos(2*pi*x)).*(cos(2*pi*y));
ue=uex(x,y)
error = (u1-ue);

%Plot error
figure, set(gcf,'DefaultAxesFontSize',8,'PaperPosition', [0 0 3.5 3.5]),  
mesh(x,y,error), colormap([0 0 0]),xlabel('x'),ylabel('y'), 
zlabel('Error'), title(strcat('Error, h=',num2str(h))); 

function [u1,x,y] = fd2poissondct(p,a,b,m)

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

% You need to deal with the zero eigenvalue here, by setting the
% appropriate value of uhat to 0.

% Computation of u = (S^(-1)*uhat)*S
u = dct(idct(uhat,1),2);

 % Append on to u the boundary values from the Neumann BC


% You should just return u, not u1.
u1 = u;

end

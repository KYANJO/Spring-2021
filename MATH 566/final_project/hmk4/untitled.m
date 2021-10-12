% Solves the Poisson equation with f(x,y) = -8*(pi^2)*cos(2*pi*x)cos(2*pi*y)


m=(2^6)-1;
a=0;b=1;
h=1/(m+1);
idx = 2:m+1;
idy = 2:m+1;
iidx=3:m+1;

[x,y] = meshgrid(a:h:b);

p=@(x,y) -8*(pi^2)*(cos(2*pi*x)).*(cos(2*pi*y));

fr = feval(p,x(idy,idx),y(idy,idx));

fhat=idct(dct(fr,1),2);

u2=zeros(m);
denom = [u2;bsxfun(@plus,cos(pi*(iidx-1)./(m+1)).',cos(pi*(iidx-1)./(m+1)))-2];

uhat = h^2/2*(fhat./denom);

u = dst(idst(uhat,1),2);

j=zeros(1,m+2);
k=zeros(m,1);
kkk=[k,u,k];
u1=[j;kkk;j];

uex=@(x,y) (cos(2*pi*x)).*(cos(2*pi*y));
ue=uex(x,y)
error = (u1-ue);
figure, set(gcf,'DefaultAxesFontSize',8,'PaperPosition', [0 0 3.5 3.5]),  
mesh(x,y,u1), colormap([0 0 0]),xlabel('x'),ylabel('y'), 
zlabel('Error'), title(strcat('Error, h=',num2str(h))); 
mesh(u1)
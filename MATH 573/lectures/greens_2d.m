function Aout = greens_2d(N)

if nargin == 0
    N = 32;
end

ic = 3;
switch ic
    case 1
        uexact = @(x,y) x + y;
        f = @(x,y) zeros(size(x));
        
    case 2
        a = 1;
        b = 1;
        uexact = @(x,y) a*(x-0.5).^2 + b*(y-0.5).^2;
        f = @(x,y) 2*(a+b)*ones(size(x));
    
    case 3
        uexact = @(x,y) sin(2*pi*x).*sin(2*pi*y);
        f = @(x,y) -8*pi^2*uexact(x,y);
end
        

ax = 0;
bx = 1;
ay = 0;
by = 1;

h = (bx-ax)/N;
xe = linspace(ax,bx,N+1);
ye = linspace(ay,by,N+1);
[xem,yem] = meshgrid(xe,ye);

% -------------------------------------------------
% Volume integral
% -------------------------------------------------
fprintf('Computing volume integral\n');
qV = zeros(N+1,N+1);
for i = 1:N+1
    for j = 1:N+1
        I = i;
        J = j;
        Px = xem(i,j);
        Py = yem(i,j);
        qV(i,j) = volume_integral(h,I,J,xem,yem,f,Px,Py);
    end
end


% -------------------------------------------------
% Solve for density mu around boundary
% -------------------------------------------------

fprintf('Computing density integral\n');

% Get arrays of (x,y) values around the bounary. 
M = 4*N;
Vx = [xe(1:end-1)'; bx*ones(N,1); fliplr(xe(2:end))'; zeros(N,1)];
Vy = [zeros(N,1); ye(1:end-1)'; by*ones(N,1); fliplr(ye(2:end))'];

% Create matrix needed to get density
A = zeros(M,M);
E = eye(M);
for k = 1:M
    Ej = E(:,k);
    for kk = 1:M
        K = kk;
        A(kk,k) = matvec_boundary(N,K,Vx,Vy,h,Ej);
    end
end

% Solve for density
% Subtract off the volume solution at the boundary. 
Vg = [qV(1:end-1,1); qV(end,1:end-1)'; ...
   flipud(qV(2:end,end)); fliplr(qV(1,2:end))'];

g1 = uexact(Vx,Vy);
g = g1 - Vg;
mu = A\g;

% -------------------------------------------------
% Get the harmonic function needed to impose boundary
% conditions. 
% -------------------------------------------------
qh = zeros(size(xem));

fprintf('Computing harmonic function\n');

% Compute boundary integral using mu found above
for i = 1:N+1
    for j = 1:N+1
        Px = xem(i,j);
        Py = yem(i,j);
        qh(i,j) = matvec_volume(N,Vx,Vy,h,Px,Py,mu);
    end
end

% Solution : Volume integral + harmonic function
q = qh + qV;


% -------------------------------------------------
% Plots
% -------------------------------------------------

% Plot the harmonic function
figure(1)
clf;


surf(xem,yem,qh);

title('Harmonic function','fontsize',18);
xlabel('x','fontsize',16);
ylabel('y','fontsize',16);
daspect([1,1,1]);


% Plot the volume integral
figure(2)
clf;

surf(xem,yem,qV);

title('Volume integral','fontsize',18);
xlabel('x','fontsize',16);
ylabel('y','fontsize',16);
daspect([1,1,1]);
axis([0,1,0,1])


% Plot the solution
figure(3)
clf;

surf(xem,yem,q);

title('Solution','fontsize',18);
xlabel('x','fontsize',16);
ylabel('y','fontsize',16);
daspect([1,1,8]);
axis([0,1,0,1])


% Plot the exact solution
figure(4)
clf;

ue = uexact(xem,yem);
surf(xem,yem,ue);

title('Exact solution','fontsize',18);
xlabel('x','fontsize',16);
ylabel('y','fontsize',16);
daspect([1,1,1]);
axis([0,1,0,1])

shg;



figure(5)
clf;

err = abs(q - ue);
surf(xem,yem,err);

title('Error','fontsize',18);
xlabel('X','fontsize',16);
ylabel('Y','fontsize',16);

h2 = h*h;
nerr(1) = sum(err(:))*h2;
nerr(2) = sqrt(sum(err(:).^2)*h2);
nerr(3) = max(err(:));
fprintf('%8d %12.4e %12.4e %12.4e %12.4f\n', ...
    N,nerr(1), nerr(2), nerr(3), cond(A));       


if nargout > 0
    Aout = A;
end

end

function Qmu = matvec_boundary(N,K,Vx,Vy,h,mu)

Px = Vx(K);
Py = Vy(K);
Rv = sqrt((Vx - Px).^2 + (Vy - Py).^2);

Rv(K) = 1e-12;

dGdn = zeros(size(Rv));
for k = 1:4
    idx = (k-1)*N + (1:N);
    if k == 1
        gradr = -(Py - Vy);
    elseif k == 2
        gradr = (Px - Vx);
    elseif k == 3
        gradr = (Py - Vy);
    else
        gradr = -(Px - Vx);
    end            
    dGdn(idx) = gradr(idx)./Rv(idx).^2;
end

Fv = -mu.*dGdn/(2*pi);
Fv(K) = 0;

Qmu = mu(K)/2 + sum(Fv)*h;


end

function Qmu = matvec_volume(N,Vx,Vy,h,Px,Py,mu)

Rv = sqrt((Vx - Px).^2 + (Vy - Py).^2);

K = find(Rv == 0);
Rv(K) = 1e-12;

dGdn = zeros(size(Rv));
for k = 1:4
    idx = (k-1)*N+(1:N);
    if k == 1
        gradr = -(Py - Vy);
    elseif k == 2
        gradr = (Px - Vx);
    elseif k == 3
        gradr = (Py - Vy);
    else
        gradr = -(Px - Vx);
    end
        
    dGdn(idx) = gradr(idx)./Rv(idx).^2;
end

Fv = -mu.*dGdn/(2*pi);

if isempty(K)
    Qmu = sum(Fv)*h;
else
    Fv(K) = 0;
    Qmu =  mu(K)/2 + sum(Fv)*h;
end

end

function fint = volume_integral(h,I,J,xem,yem,f,Px,Py)

rem = sqrt((xem - Px).^2 + (yem - Py).^2);

% Avoid the singularity
rem(I,J) = 1e-12;

F = f(xem,yem).*log(rem)/(2*pi);
F_int = sum(sum(F(2:end-1,2:end-1)));

F_edge = 0.5*(sum(F(2:end-1,1))   + sum(F(end,2:end-1)) + ...
              sum(F(2:end-1,end)) + sum(F(1,  2:end-1)));
   
F_corners = 0.25*(F(1,1) + F(end,1) + F(end,end) + F(1,end));

fint = (F_int + F_edge + F_corners)*h*h;

end

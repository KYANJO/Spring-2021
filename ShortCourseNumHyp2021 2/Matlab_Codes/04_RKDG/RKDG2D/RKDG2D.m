%************************************************************************** 
%% Unstrutured RK-DG method for 2D hyperbolic conservation laws
% Q_t + \nabla \cdot F = 0 with F = (A1,A2)*Q
clear all
close all
clc

global U0 V0 c
global nVar nElem N Neighbor NeighborEdge NormalVector Area EdgeLength
global nDOF M Kxi Keta Fm Fp x y tri

%% Model related parameters
c = 1;        
nVar = 4;
U0 = 0; % Advection speed in x direction
V0 = 0; % Advection speed in y direction

% Runge-Kutta scheme
rktype = 4;
% Polynomial approximation degree
N= 3;

%% Define computational domain
xL = -1;                % Left boundary
xR = 1;                 % Right boundary
yL = -1;                % Bottom boundary
yR = 1;                 % Top boundary
IMAX = 10;              % Number of cells in x-direction (boundary)
JMAX = 10;              % Number of cells in y-direction (boundary)
dx = (xR-xL)/IMAX;      % Mesh spacing on the boundary
dy = (yR-yL)/JMAX;
nNode = 0;              % Counter for the total number of nodes

%% Define data for time discretization 
CFL = 0.45;             % Courant-Friedrichs-Lewy number
time = 0;               % Initial time
tend = 0.5;             % End time
NMAX = 1000;            % Maximum number of time steps

% Edge definition for each triangle
EdgeDef = [1 2;             % edge 1 of triangle i
    2 3;                    % edge 2
    3 1];                   % edge 3
% Initialize the necessary data for DG
InitDG;

%% Create mesh
% Point cloud
for i = 1:IMAX+1
    for j= 1:JMAX+1
        nNode = nNode +1;
        x(nNode) = xL + (i-1)*dx;   % x coordinate of the node
        y(nNode) = yL + (j-1)*dy;   % y coordinate of the node
    end
end
tri = delaunay(x,y);                % Compute the delalunay triangulation
[nElem,tmp] = size(tri);            % Compute the number of elements (nElem)
triplot(tri,x,y)                    % Plot the mesh

% Compute the elements attached to each node (connectivity of 
% the dual mesh, Voronoi teselation)
nVertexElement = zeros(nNode,1); % Number of elements attached to a vertex 
VertexElement = zeros(nNode,10); % List of elements attached to avertex
for i=1:nElem   % Loop on elements
    for k=1:3       % Loop over the 3 nodes of a triangle
        iNode = tri(i,k);        % Recover the global index of the node
        nVertexElement(iNode) = nVertexElement(iNode) + 1;
        VertexElement(iNode,nVertexElement(iNode)) = i;
    end
end

% Compute the neighbors of each triangle
Neighbor = zeros(nElem,3);      % each triangle can have 3 face neighbors
NeighborEdge = zeros(nElem,3);  % local edge number inside the neighbor (not used in FV)
% Search neighbors
for i=1:nElem
    for iEdge=1:3
        iNode1= tri(i,EdgeDef(iEdge,1)); % global node number
        iNode2= tri(i,EdgeDef(iEdge,2));
        for k=1:nVertexElement(iNode1)
            j = VertexElement(iNode1,k);
            if (i==j) % element i cannot be its own enighbor
                continue
            end
            for jEdge = 1:3
                jNode1  = tri(j,EdgeDef(jEdge,1));
                jNode2  = tri(j,EdgeDef(jEdge,2));
                if (iNode1 == jNode2 && iNode2 == jNode1)
                    Neighbor(i,iEdge) = j;
                    NeighborEdge(i,iEdge) = jEdge; % not used in FV
                    break
                end
            end
            if (Neighbor(i,iEdge)>0)
                break % neighbor has been found, exit loop
            end
        end
    end
end
% Compute the normal vectors, areas, edge lengths
z = [0;0;1]; %vector e_z
for i=1:nElem
    % Normal vectors and edge lengths
    check = zeros(2,1);  % check mesh for consistency
    for iEdge = 1:3
        iNode1= tri(i,EdgeDef(iEdge,1)); % global node number
        iNode2= tri(i,EdgeDef(iEdge,2));
        v(:,iEdge) = [x(iNode2)-x(iNode1);
            y(iNode2)-y(iNode1);
            0];
        n3d = cross( v(:,iEdge),z); % outward normal vector
        NormalVector(:,iEdge,i) = n3d(1:2)/sqrt(sum(n3d.^2)); % Unit normal vector
        EdgeLength(iEdge,i) = sqrt(sum(v(:,iEdge).^2));
        check = check + NormalVector(:,iEdge,i)*EdgeLength(iEdge,i);
    end
    if (sqrt(sum(check.^2))>1e-12)
        disp('Mesh is not consistent. Please debug!')
        return
    end
    % Area
    temp = cross(v(:,1),-v(:,3));
    Area(i) = temp(3)/2;
    % Barycenter of each triangle
    xb(i) = 1/3*sum(x(tri(i,:)));
    yb(i) = 1/3*sum(y(tri(i,:)));    
    % Incircle diameter
    Incircle(i) = 4*Area(i)/sum(EdgeLength(:,i));
end
disp('Mesh seems to  be consistent')

% Create subtriangulation for plotting
nNodePlot = nDOF*nElem;
nElemPlot = nSub*nElem;
xp = zeros(nNodePlot,1);
yp = zeros(nNodePlot,1);
trip = zeros(nElemPlot,3);
countE = 0; % element counter
countN = 0; % node counter
for iElem=1:nElem
    % create subtriangulation for the plotting of the high order DG
    % polynomials
    offset = (iElem-1)*nDOF;
    for iSub=1:nSub
        countE = countE + 1;
        trip(countE,:) = offset + SubTri(iSub,:);
    end
    % create the subnodes
    X = x(tri(iElem,:));
    Y = y(tri(iElem,:));
    for iSubNode=1:nDOF
        countN = countN + 1;
        xi  = xisub(iSubNode,1);
        eta = xisub(iSubNode,2);
        xx = X(1) + xi *(X(2)-X(1)) + eta*(X(3)-X(1));
        yy = Y(1) + xi *(Y(2)-Y(1)) + eta*(Y(3)-Y(1));
        xp(countN) = xx;
        yp(countN) = yy;
    end
end
disp('Subgrid triangulation created')

%% Set initial condition L2 projection
uhat = zeros(nDOF,nVar,nElem);
for i=1:nElem
    X = x(tri(i,:));
    Y = y(tri(i,:));
    for iGP = 1:nGP
        % Coordinate Gaussian quadrature point in reference coordinates
        xi = xiGP(iGP,1);
        eta = xiGP(iGP,2);
        % Physical coordinates
        xGP = X(1) + xi*(X(2)-X(1)) +eta*(X(3)-X(1));
        yGP = Y(1) + xi*(Y(2)-Y(1)) +eta*(Y(3)-Y(1));
        % Define IC
        Q = [1;0;0;0];
        Q(4) = Q(4) + exp(-0.5*(xGP^2+yGP^2)/0.2^2);
        % Compute basis function in the quadrature point
        [phi,phi_xi,phi_eta] = BaseFunc(xi,eta);
        for k=1:nDOF
            uhat(k,:,i) = uhat(k,:,i) + wGP(iGP)*phi(k)*Q(:)';
        end
    end
    % Multiply by the inverse of the mass matrix
    % Our basis is orthogonal
    for k=1:nDOF
        uhat(k,:,i) = uhat(k,:,i)/M(k,k);
    end
end
% Plot IC
Qplot = zeros(nVar,nNodePlot);
count = 0;
for i=1:nElem
    for iSubNode = 1:nDOF
        count = count + 1;
        xi = xisub(iSubNode,1);
        eta = xisub(iSubNode,2);
        [phi,phi_xi,phi_eta] = BaseFunc(xi,eta);
        Qplot(:,count) = (phi(:)'*uhat(:,:,i))';
    end
end
s = trisurf(trip,xp,yp,Qplot(4,:));
set(s,'EdgeColor',[0,0,0])
set(s,'FaceColor','interp')


tic
for n=1:NMAX
    % Compute the maximum abs eigenvalues
    amax = 0;
    for i=1:nElem
        amax = max(amax, max(abs(lambda(uhat(1,:,i),[0,0]))) );
    end
    dt = CFL*min(Incircle)/(amax*(2*N+1));
    if (time +dt >tend)
        dt = tend-time;
    end
    if (time >= tend)
        break
    end
    
    % RK scheme
    switch(rktype)
        case(1)
            % Forward Euler 
            k1 = LhDG(uhat);
            uhatnew = uhat + dt*k1; 
        case(2)
            % RK2 
            k1 = LhDG(uhat);
            k2 = LhDG(uhat+dt/2*k1); 
            uhatnew = uhat + dt*k2; 
        case(4)
            % RK4 
            k1 = LhDG(uhat);
            k2 = LhDG(uhat+dt/2*k1);
            k3 = LhDG(uhat+dt/2*k2);
            k4 = LhDG(uhat+  dt*k3);
            uhatnew = uhat + dt/6*(k1+2*k2+2*k3+k4);
    end
    % Update time and overwrite solution
    time = time+dt;
    uhat = uhatnew;
    
    % Plot solution
    Qplot = zeros(nVar,nNodePlot);
    count = 0;
    for i=1:nElem
        for iSubNode = 1:nDOF
            count = count + 1;
            xi = xisub(iSubNode,1);
            eta = xisub(iSubNode,2);
            [phi,phi_xi,phi_eta] = BaseFunc(xi,eta);
            Qplot(:,count) = (phi(:)'*uhat(:,:,i))';
        end
    end
    s = trisurf(trip,xp,yp,Qplot(4,:));
    set(s,'EdgeColor',[0,0,0])
    set(s,'FaceColor','interp')
    title(sprintf('Time = %f',time));
    drawnow
end
toc
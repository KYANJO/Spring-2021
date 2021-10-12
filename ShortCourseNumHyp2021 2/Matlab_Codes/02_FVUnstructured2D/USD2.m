%% Unstrutured FV method for 2D hyperbolic conservation laws
% Q_t + \nabla \cdot F = 0
clear all
close all
clc

global gamma grav nVar eqntype

%% Model related parameters
eqntype = 2;
switch(eqntype)
    case(1)% SW
        grav = 9.81;
        nVar = 3;
        ivarplot = 1;        
    case(2) % Euler
        gamma = 1.4;        
        nVar = 4;
        ivarplot = 4;
end

%% Define computational domain
xL = -1;                % Left boundary
xR = 1;                 % Right boundary
yL = -1;                % Bottom boundary
yR = 1;                 % Top boundary
IMAX = 20;              % Number of cells in x-direction (boundary)
JMAX = 20;              % Number of cells in y-direction (boundary)
dx = (xR-xL)/IMAX;      % Mesh spacing on the boundary
dy = (yR-yL)/JMAX;
nNode = 0;              % Counter for the total number of nodes

%% Define data for time discretization 
CFL = 0.45;             % Courant-Friedrichs-Lewy number
time = 0;               % Initial time
tend = 1.0;             % End time
NMAX = 1000;            % Maximum number of time steps

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
% Edge definition for each triangle
EdgeDef = [1 2;             % edge 1 of triangle i
    2 3;                    % edge 2
    3 1];                   % edge 3
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

%% Define initial condition
for i=1:nElem
    switch(eqntype)
        case(1)
            % Shallow water
            Q(:,i) = [ 1; 0; 0];
            Q(1,i) = Q(1,i) + 1*exp(-0.5*( xb(i)^2+yb(i)^2 )/0.2^2);
        case(2)
            % Euler
            Q(:,i) = [ 1; 0; 0; 1/(gamma-1) ];
            Q(4,i) = Q(4,i) + 1*exp(-0.5*( xb(i)^2+yb(i)^2 )/0.2^2);
    end
end

% Interpolate data from cell centers to the nodes
NodeArea = zeros(nNode,1);
QNode    = zeros(nVar,nNode);
for i=1:nElem
    for k=1:3
        iNode = tri(i,k);
        NodeArea(iNode) = NodeArea(iNode) + Area(i);
        QNode(:,iNode)  = QNode(:,iNode)  + Area(i)*Q(:,i);
    end
end
% area-weighted average
for iNode=1:nNode
        QNode(:,iNode) = QNode(:,iNode) / NodeArea(iNode);
end

% Plot IC
s = trisurf(tri,x,y,QNode(ivarplot,:));
set(s,'EdgeColor',[0,0,0])  % black edges
set(s,'FaceColor','interp')

%% FV solver
for n=1:NMAX % Time loop
    % Compute maximum eigenvalue
    amax = 0;    
    for i=1:nElem
        amax = max(amax, max(abs(lambda(Q(:,i),[0;0]))) ); 
    end 
    % Compute time step
    dt = CFL*min(Incircle)/amax;
    if(time+dt>tend)
        dt = tend-time;
    end
    if(time>=tend)
        break
    end
    Qnew = Q;
    for i=1:nElem
        for iEdge=1:3
            j = Neighbor(i,iEdge); % find the index of the element on the other side
            nij = NormalVector(:,iEdge,i); % normal vector to the face
            if (j==0)
                % we are in the boundary: reflective wall
                Qi = Q(:,i);
                Qj = Qi;
                Qj(2:3) = Qj(2:3) -2 * (nij(1)*Qi(2)+nij(2)*Qi(3))*nij;
            else
                % Inside thedomain
                Qi = Q(:,i);
                Qj = Q(:,j);
            end
            smax = max( max(abs(lambda(Qi,nij))), max(abs(lambda(Qj,nij))) );
            % Numerical flux in normal direction (Rusanov/local Lax-Fiedrichs flux)
            Fn = 0.5*(F(Qi)+F(Qj))*nij - 0.5 *smax*(Qj-Qi);
            Qnew(:,i) = Qnew(:,i) - dt/Area(i)*EdgeLength(iEdge,i)*Fn;
        end
    end
    
    % Update time and overwrite solution
    time = time +dt;
    Q = Qnew;
    
    % Interpolate from cell centers to nodes
    QNode    = zeros(nVar,nNode);
    for i=1:nElem
        for k=1:3
            iNode = tri(i,k);
            QNode(:,iNode)  = QNode(:,iNode)  + Area(i)*Q(:,i)/ NodeArea(iNode);
        end
    end
    % Plot solution 
    s = trisurf(tri,x,y,QNode(ivarplot,:));
    set(s,'EdgeColor',[0,0,0])  % black edges
    set(s,'FaceColor','interp')
    title(sprintf('Current time = %f',time))
%     axis([xL xR yL yR 0 2])
    drawnow
end
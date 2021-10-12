% Initialize the RKDG scheme 
% 
nDOF = (N+1)*(N+2)/2;                   % compute the number of degrees of freedom on triangles 
[nGP, wGP, xiGP] = TriQuadPoints;       % provide the quadrature points on triangles (valid up to N=3). 
[chi1D,w1D]=gauleg(N+1);                % compute the 1D Gauss-Legendre quadrature points and weights 
%
% elementary volume integrals 
% 
% compute the mass matrix and the spatial stiffness matrices using Gaussian
% quadrature (exact for polnomials up to degree N=3) 
M     = zeros(nDOF,nDOF);    % mass matrix 
Kxi   = zeros(nDOF,nDOF);    % stiffness matrix containing the xi derivative 
Keta  = zeros(nDOF,nDOF);    % stiffness matrix containing the eta derivative 
ddxi  = zeros(nDOF,nDOF);    % spatial xi derivatives 
ddeta = zeros(nDOF,nDOF);    % spatial eta derivatives 
% Loop over the 2D quadrature points on the triangle 
for iGP=1:nGP
    xi = xiGP(iGP,1); 
    eta = xiGP(iGP,2); 
    [phi, phi_xi, phi_eta] = BaseFunc(xi,eta);   
    for k=1:nDOF  
        for m=1:nDOF  
            M(k,m)    = M(k,m)    + wGP(iGP)*phi(k)*phi(m); 
            Kxi(k,m)  = Kxi(k,m)  + wGP(iGP)*phi_xi(k)*phi(m); 
            Keta(k,m) = Keta(k,m) + wGP(iGP)*phi_eta(k)*phi(m); 
        end
    end
end
% matrices which provide the polyomials of the spatial derivatives in xi
% and eta direction 
ddxi  = inv(M)*Kxi'; 
ddeta = inv(M)*Keta'; 
% elementary boundary integrals 
% 
% compute the matrices that integrate the numerical flux on the boundary 
% quadrature (exact for polnomials up to degree N=3)  
Fm    = zeros(nDOF,nDOF,3);    % flux matrices from within the element
Fp    = zeros(nDOF,nDOF,3,3);  % flux matrices from outside the element (all possible configurations need to be stored!) 
% The vertices of the unit element 
UnitElem = [0,0;
            1,0;
            0,1]; 
X1=zeros(2,1);
X2=zeros(2,1);
for iSide=1:3    
    pNode1 = EdgeDef(iSide,1);      % get node number on the unit element of the first node on the edge 
    pNode2 = EdgeDef(iSide,2);      % get node number on the unit element of the second node on the edge 
    % coordinates of the two nodes 
    X1(1) = UnitElem(pNode1,1); 
    X1(2) = UnitElem(pNode1,2); 
    X2(1) = UnitElem(pNode2,1); 
    X2(2) = UnitElem(pNode2,2);    
    % loop over the 1D Gauss-Legendre points along the edge 
    for iGP=1:N+1
        chi = chi1D(iGP); 
        xi  = X1(1) + chi*(X2(1)-X1(1)); 
        eta = X1(2) + chi*(X2(2)-X1(2)); 
        [phi, phi_xi, phi_eta] = BaseFunc(xi,eta); 
        for k=1:nDOF
            for m=1:nDOF
                Fm(k,m,iSide) = Fm(k,m,iSide) + w1D(iGP)*phi(k)*phi(m); 
            end
        end
    end
    % for the Fp matrix, loop over all possible configurations in the
    % neighbor element
    for jSide=1:3
        qNode1 = EdgeDef(jSide,1);
        qNode2 = EdgeDef(jSide,2); 
        % coordinates of the two nodes 
        Y1(1) = UnitElem(qNode1,1); 
        Y1(2) = UnitElem(qNode1,2); 
        Y2(1) = UnitElem(qNode2,1); 
        Y2(2) = UnitElem(qNode2,2);    
        % loop over the 1D Gauss-Legendre points along the edge 
        for iGP=1:N+1
            % in the element itself 
            chi = chi1D(iGP); 
            xi  = X1(1) + chi*(X2(1)-X1(1)); 
            eta = X1(2) + chi*(X2(2)-X1(2)); 
            [phi, phi_xi, phi_eta] = BaseFunc(xi,eta); 
            % in the neighbor. counterclockwise orderning in the element =>
            % clockwise in the neighbor !!
            zeta = 1-chi; 
            xin  = Y1(1) + zeta*(Y2(1)-Y1(1)); 
            etan = Y1(2) + zeta*(Y2(2)-Y1(2)); 
            [phin, phi_xi, phi_eta] = BaseFunc(xin,etan); 
            for k=1:nDOF
                for m=1:nDOF
                    Fp(k,m,iSide,jSide) = Fp(k,m,iSide,jSide) + w1D(iGP)*phi(k)*phin(m); 
                end
            end
        end        
    end
end
%
% for the data output create a subtriangulation 
% 
xisub = zeros(nDOF,2);      % coordinates of the subtriangles 
nSub  = N^2;                % number of subtriangles 
SubTri = zeros(nSub,3);     % sub-connectivity 
count = 0; 
idx = zeros(N+1,N+1); 
for j=0:N
    for i=0:N-j 
        count = count + 1; 
        idx(i+1,j+1) = count; 
        xisub(count,:) = [i/N, j/N]; 
    end
end
count = 0; 
% upward pointing triangles 
for j=1:N+1
    for i=1:N+1-j
        count = count + 1;
        iNode1 = idx(i,j); 
        iNode2 = idx(i+1,j); 
        iNode3 = idx(i,j+1); 
        SubTri(count,:) = [iNode1,iNode2,iNode3]; 
    end
end
% downward pointing triangles                 
for j=1:N
    for i=2:N+1-j
        count = count + 1;
        iNode1 = idx(i,j); 
        iNode2 = idx(i,j+1); 
        iNode3 = idx(i-1,j+1); 
        SubTri(count,:) = [iNode1,iNode2,iNode3]; 
    end
end
    
% % unit test for the compatibility of the stiffness and flux matrices 
% fhat = zeros(nDOF,1); 
% ghat = zeros(nDOF,1); 
% fhat(1) = 1; 
% ghat(1) = 1; 
% test  = -(Kxi*fhat + Keta*ghat); 
% testj = zeros(nDOF,3); 
% for jEdge=1:3
%     testj(:,jEdge) = -(Kxi*fhat + Keta*ghat); 
% end
% % 
% for iEdge=1:3
%     pNode1 = EdgeDef(iEdge,1);      % get node number on the unit element of the first node on the edge 
%     pNode2 = EdgeDef(iEdge,2);      % get node number on the unit element of the second node on the edge 
%     % coordinates of the two nodes 
%     X1(1) = UnitElem(pNode1,1); 
%     X1(2) = UnitElem(pNode1,2); 
%     X2(1) = UnitElem(pNode2,1); 
%     X2(2) = UnitElem(pNode2,2);    
%     vv = X2-X1;                     % tangent vector 
%     elen = sqrt(vv(1)^2+vv(2)^2);  
%     nij = [ vv(2); -vv(1) ]/elen;   % unit normal vector 
%     test = test + Fm(:,:,iEdge)*( fhat(:)*nij(1) + ghat(:)*nij(2) )*elen; 
%     for jEdge=1:3
%         testj(:,jEdge) = testj(:,jEdge) + Fp(:,:,iEdge,jEdge)*( fhat(:)*nij(1) + ghat(:)*nij(2) )*elen; 
%     end
% end
% 
% if(norm(test)>1e-13)
%     disp('Warning. Consistency test was not successful.')
%     test
%     pause
% end
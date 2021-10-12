% Spatial discretization operator of a quadrature-free (!!) Runge-Kutta DG scheme 
% extension to nonlinear PDE systems requires the use of a nodal basis and
% therefore the elementary mass matrix will in general not be orthogonal
% any more, but this is no problem for our linear case. 
function duhatdt = LhDG(uhat) 
global nElem nVar Neighbor NeighborEdge NormalVector Area EdgeLength x y tri 
global N nDOF M Kxi Keta Fm Fp   

% Compute flux
fhat = zeros(nDOF,nVar,nElem);
ghat = zeros(nDOF,nVar,nElem);
for i=1:nElem
    for k=1:nDOF
        [f g] = flux(uhat(k,:,i)); % Physical flux
        fhat(k,:,i) = f';
        ghat(k,:,i) = g';
    end
end

duhatdt = zeros(nDOF,nVar,nElem);
for i=1:nElem
    %% Compute the Jacobian matrix of the mapping 
    X = x(tri(i,:));
    Y = y(tri(i,:));
    J = [X(2)-X(1), X(3)-X(1);
        Y(2)-Y(1), Y(3)-Y(1)];
    detJ = J(1,1)*J(2,2)-J(1,2)*J(2,1);
    iJ = [J(2,2), -J(1,2);
        -J(2,1), J(1,1)];
    iJ = iJ/detJ;
    
    % Compute the integral on the volume
    duhatdt(:,:,i) = (Kxi*(iJ(1,1)*fhat(:,:,i) + iJ(1,2)*ghat(:,:,i)) ...
        + Keta*( iJ(2,1)*fhat(:,:,i) + iJ(2,2)*ghat(:,:,i) ));
    
    % Compute flux in the boundary
    for iEdge = 1:3
        j = Neighbor(i,iEdge);
        jEdge = NeighborEdge(i,iEdge);
        nij = NormalVector(:,iEdge,i);
        if (j==0) % Boundary edge
            Fn = ( Fm(:,:,iEdge)*(fhat(:,:,i)*nij(1)+ ghat(:,:,i)*nij(2)) );
        else
            % Inside of the domain
            smax = max( max(abs(lambda(uhat(1,:,i),nij))), max(abs(lambda(uhat(1,:,j),nij))) );
            % Projection of Rusanov flux in normal direction
            Fn = 0.5*(Fp(:,:,iEdge,jEdge)*(fhat(:,:,j)*nij(1)+ghat(:,:,j)*nij(2))...
                + Fm(:,:,iEdge)*(fhat(:,:,i)*nij(1)+ghat(:,:,i)*nij(2))) ...
                - 0.5*smax*( Fp(:,:,iEdge,jEdge)*uhat(:,:,j) - Fm(:,:,iEdge)*uhat(:,:,i) );
        end
        duhatdt(:,:,i) = duhatdt(:,:,i) - 1/detJ*EdgeLength(iEdge,i)*Fn;
    end
    % "Multiply" by the inverse of the elementary mass matrix
    for k=1:nDOF
        duhatdt(k,:,i) = duhatdt(k,:,i)/M(k,k);
    end
end
end
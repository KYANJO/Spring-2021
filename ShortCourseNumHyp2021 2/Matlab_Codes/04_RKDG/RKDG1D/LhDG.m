function dqdt = LhDG(uhat)
global nGP wGP xiGP; 
global iMe Me phiL phiR phiGP phi_xiGP; 
global N dx IMAX; 
% Initialize
dqdt = zeros(N+1,IMAX);

%% Compute boundary extrapolated values
for i=1:IMAX
    wR(i) = phiR(:)'*uhat(:,i);% boundary-extrapolated value of the DG scheme for the right interface 
    wL(i) = phiL(:)'*uhat(:,i);% boundary-extrapolated value of the DG scheme for the left interface  
end

%% Compute numerical fluxes and volume integral
for i=2:IMAX-1
    % Fluxes
    smax = max( abs( a(wL(i+1)) ), abs( a(wR (i)) ) );
    fp = 0.5*( f(wR(:,i)) + f(wL(:,i+1)) ) - 0.5*smax*(wL(:,i+1) - wR(:,i));
    smax = max( abs( a(wL(i)) ), abs( a(wR (i-1)) ) );
    fm = 0.5*( f(wR(:,i-1)) + f(wL(:,i)) ) - 0.5*smax*(wL(:,i) - wR(:,i-1));
    for k=1:N+1
        dqdt(k,i) = dqdt(k,i) + -( phiR(k)*fp -phiL(k)*fm );
    end
    
    % Volume integrals
    for j=1:nGP
        uh = phiGP(:,j)'*uhat(:,i);
        fh = f(uh);
        for k=1:N+1
            dqdt(k,i) = dqdt(k,i) + wGP(j)*phi_xiGP(k,j)*fh;
        end
    end
    
    dqdt(:,i) = 1/dx*iMe*dqdt(:,i);
end
end
% Routine to initialize Gaussian points, Galerkin matrices...
%************************************************************
%% Define Gaussian quadrature points
nGP = N+1;
[xiGP,wGP] = gauleg(nGP);

%% Initialize space-time DG matrices
M = zeros(nDOF,nDOF);       % theta*theta
Kxi = zeros(nDOF,nDOF);     % theta*theta_xi
Ktau = zeros(nDOF,nDOF);   % theta*theta_tau
F0 = zeros(nDOF,nDOFs);     % space integral t=0
F1 = zeros(nDOF,nDOF);       % space integral t=1
%
Fp = zeros(1,nDOF);              % Right boundary time integral
Fm = zeros(1,nDOF);             % Left boundary time integral
% 
MSrc = zeros(nDOFs,nDOF);   % psi*theta

%% Prepare fine output
cnt = 1; % auxiliary counter
% Compute coordinates
for i=1:IMAX
    for j = 1:NFINE+1
        xi = (j-1)/NFINE;
        xfine(cnt) = x(i) - dx/2 +xi*dx;
        cnt = cnt+1;
    end
end
% Compute base function at fine output points
for j= 1:NFINE+1
    if (N==0)
        xi = 0.5;
    else
        xi = (j-1)/NFINE;
    end
    basefine(j,:) = SpaceBaseFunc(xi);
end
%% Compute basis functions at the Gaussian points
for j = 1:nGP
    xi = xiGP(j);
    basegp(j,:) = SpaceBaseFunc(xi);
end
%% Compute space-time integrals
disp('Computing space-time integrals...')
for itGP = 1:nGP
    for ixGP = 1:nGP
        xi  = xiGP(ixGP);
        tau = xiGP(itGP);
        weight = wGP(itGP)*wGP(ixGP);
        [theta,theta_xi,theta_tau]  = SpaceTimeBaseFunc(xi,tau);
        [psi,psi_xi]                = SpaceBaseFunc(xi);
        for k=1:nDOF
            for l=1:nDOF
                M(k,l)     = M(k,l)     + weight*theta(k)*theta(l);
                Kxi(k,l)   = Kxi(k,l)   + weight*theta(k)*theta_xi(l);
                Ktau(k,l)  = Ktau(k,l)  + weight*theta_tau(k)*theta(l);
            end
        end
        for k=1:nDOFs
            for l=1:nDOF
                MSrc(k,l)   = MSrc(k,l)   + weight*psi(k)*theta(l); 
            end
        end            
    end
end
Kx   = Kxi';
iM   = inv(M); 
dudx = iM*Kxi;  
%% Compute space integrals at t=1 and t=0
for ixGP = 1:nGP
    xi  = xiGP(ixGP);
    tau = 1;
    [theta1,theta1_xi,theta1_tau]=SpaceTimeBaseFunc(xi,tau);
    tau = 0;
    [theta0,theta0_xi,theta0_tau]=SpaceTimeBaseFunc(xi,tau);
    for k=1:nDOF
        for l=1:nDOF
            weight = wGP(ixGP);
            F1(k,l)    = F1(k,l)    + weight*theta1(k)*theta1(l);
        end
    end
    [psi,psi_xi] = SpaceBaseFunc(xi);
    for k=1:nDOF
        for l=1:nDOFs
            weight = wGP(ixGP);
            F0(k,l)    = F0(k,l)    + weight*theta0(k)*psi(l);
        end
    end
end
%% Compute boundary time integrals
for itGP = 1:nGP
    weight = wGP(itGP);
    tau = xiGP(itGP);
    xi  = 0;
    [theta0,theta0_xi,theta0_tau]   = SpaceTimeBaseFunc(xi,tau);
    xi  = 1;
    [theta1,theta1_xi,theta1_tau]   = SpaceTimeBaseFunc(xi,tau);
    for l=[1:nDOF]
        Fp(l)   = Fp(l)   + weight*1*theta1(l);
        Fm(l)   = Fm(l)   + weight*1*theta0(l);
    end
end

%% Matrices
K1  = F1 - Ktau;                     
iK1 = inv(K1);        







function lambda = PDEEigenvalues(Q)
%% Function to compute the eigenvalues of the Jacobian
% matrix of the system
global pdetype gamma g
switch(pdetype)
    case(0) % Euler
        p = (gamma-1)*( Q(3) - 0.5*Q(2)^2/Q(1) );
        c = sqrt(gamma*p/Q(1));
        vel = Q(2)/Q(1);
        lambda(1,1) = vel-c;
        lambda(2,1) = vel;
        lambda(3,1) = vel+c;
    case(1) % Shallow water equations with scalar transport 
        c = sqrt(g*Q(1));
        vel = Q(2)/Q(1);
        lambda(1,1) = vel-c;
        lambda(2,1) = vel;
        lambda(3,1) = vel+c;
    case(2) % Shallow water equations with variable bottom
        c = sqrt(g*Q(1));
        vel = Q(2)/Q(1);
        lambda(1,1) = vel-c;
        lambda(2,1) = vel;
        lambda(3,1) = vel+c;
end

end
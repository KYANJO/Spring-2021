function L = PDEEigenvalues(q)
% Function containing the eigenvalues of our system
global pdetype gamma
switch(pdetype)
    case(0) % Euler equations
        vel = q(2)/q(1);
        p = (gamma-1)*(q(3)-0.5*q(2)^2/q(1));
        c = sqrt(gamma*p/q(1));
        
        L(1,1) = vel-c;
        L(2,1) = vel;
        L(3,1) = vel+c;
end
end
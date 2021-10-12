function f = PDEFlux(Q)
% Compute the nonlinear flux as a function of the
% vector of conservative variables Q
global pdetype gamma g

switch(pdetype)
    case(0) % Euler equations
        p = (gamma-1)*( Q(3) - 0.5*Q(2)^2/Q(1) );
        f(1) = Q(2);
        f(2) = Q(2)^2/Q(1) + p;
        f(3) = Q(2)/Q(1)*( Q(3) + p );
    case(1) % Shallow water equations with scalar transport  
        f(1) = Q(2);
        f(2) = Q(2)^2/Q(1) + 0.5*g*Q(1)^2;
        f(3) = Q(2)*Q(3)/Q(1);
    case(2) % Shallow water equations with variable bottom
        f(1) = Q(2);
        f(2) =  Q(2)^2/Q(1) + 0.5*g*Q(1)^2;
        f(3) = 0;
end

end
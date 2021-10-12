function V = PDECons2Prim(Q)
% Compute primitive variables from conservative variables
global pdetype gamma
switch(pdetype)
    case(0) % Euler equations
        V(1) = Q(1);                         % Density
        V(2) = Q(2)/Q(1);                 % Velocity
        V(3) = (gamma-1)*(Q(3) - 0.5*Q(2)^2/Q(1));  % Pressure
    case(1) % Shallow water equations with scalar transport 
        V(1) = Q(1);                         
        V(2:3) = Q(2:3)/Q(1);                 
    case(2) % Shallow water equations with variable bottom
        V(1) = Q(1);                       
        V(2) = Q(2)/Q(1); 
        V(3) = Q(3);
end
end
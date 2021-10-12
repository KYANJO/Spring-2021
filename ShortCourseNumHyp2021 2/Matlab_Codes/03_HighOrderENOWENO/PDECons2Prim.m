function V = PDECons2Prim(Q)
% Compute primitive variables from given conservative variables
global pdetype gamma
switch(pdetype)
    case(0)
        V(1) = Q(1); % Density
        V(2) = Q(2)/Q(1); % Velocity
        V(3) = (gamma-1)*( Q(3) - 0.5*Q(2)^2/Q(1) ); % Pressure
end
end
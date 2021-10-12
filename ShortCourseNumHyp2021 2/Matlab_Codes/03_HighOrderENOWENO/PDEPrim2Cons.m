function Q = PDEPrim2Cons(V)
% Compute conservative variables from given primitive variables
global pdetype gamma
switch(pdetype)
    case(0)
        Q(1) = V(1); % Density
        Q(2) = V(1)*V(2); % Linear momentum
        Q(3) = V(3)/(gamma-1) + 0.5*V(1)*V(2)^2; % Total energy density
end
end
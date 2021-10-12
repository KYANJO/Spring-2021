function Q = PDEPrim2Cons(V)
% Compute the conservative variables from primitive variables
global pdetype gamma
switch(pdetype)
    case(0)
        Q(1) = V(1);                         % Density
        Q(2) = V(1)*V(2);                 % Linear momentum
        Q(3) = V(3)/(gamma-1) + 0.5*V(1)*V(2)^2;  % Total energy
    case(1)
        Q(1) = V(1); 
        Q(2:3) = V(1)*V(2:3); 
    case(2)
        Q(1) = V(1);
        Q(2) = V(1)*V(2);
        Q(3) = V(3);
end
end
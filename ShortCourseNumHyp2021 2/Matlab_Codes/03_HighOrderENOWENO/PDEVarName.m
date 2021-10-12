function Name = PDEVarName
% Function to define the name of the primitive variables f the system
global pdetype
switch(pdetype)
    case(0) % Euler equations
        Name(1,1:3) = 'rho';
        Name(2,1) = 'u';
        Name(3,1) = 'p';
end
end
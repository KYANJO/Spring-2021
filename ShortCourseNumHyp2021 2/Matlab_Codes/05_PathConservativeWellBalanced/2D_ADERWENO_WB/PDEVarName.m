function varname=PDEVarName
global pdetype
% names of the primitive variables 
switch(pdetype)
    case(0) % Euler equations 
        varname(1,1:3) = 'rho';
        varname(2,1) = 'u';
        varname(3,1) = 'p';
    case(1) % Shallow water equations with scalar transport 
        varname(1,1) = 'h';
        varname(2,1) = 'u';
        varname(3,1:3) = 'psi';
    case(2) % Shallow water equations with variable bottom 
        varname(1,1) = 'h';
        varname(2,1) = 'u';
        varname(3,1) = 'b';
end
end
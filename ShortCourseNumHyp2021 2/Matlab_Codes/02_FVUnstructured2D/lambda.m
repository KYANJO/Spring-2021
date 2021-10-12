%% Eigenvalues of the PDE 
function L=lambda(Q,n) 
global gamma grav eqntype
u = Q(2)/Q(1); 
v = Q(3)/Q(1);
switch (eqntype)
    case(1)
         c = sqrt(grav*Q(1));                        % SWE 
    case(2)
        p = (gamma-1)*( Q(4) - 0.5*Q(1)*(u^2+v^2) ); % Euler 
        c = sqrt(gamma*p/Q(1));                      
end

if( abs(n(1))+abs(n(2)) == 0 )
    un = sqrt(u^2+v^2);
else
    un = u*n(1) + v*n(2); 
end

% Eigenvalues in normal direction   
switch (eqntype)
    case(1) %SWE
        L=[un-c, un, un+c];
    case(2)% Euler
        L=[un-c, un, un, un+c];  
end

end
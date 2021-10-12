function flux=f(q)
%Flux function

% Linear scalar advection equation 
%flux = +q; 

% Burgers equation 
flux = q^2/2; 

% cubic flux 
%flux = q^3/3; 
end
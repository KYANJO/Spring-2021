function dfdq=a(q)
% Function with the jacobian of the flux

% Linear scalar advection equation 
%dfdq = +1; 

% Burgers equation 
dfdq = q; 

% cubic flux 
%dfdq = q^2; 
end
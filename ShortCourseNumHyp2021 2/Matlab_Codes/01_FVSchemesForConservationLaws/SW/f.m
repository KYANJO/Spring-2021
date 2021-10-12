% Function of the physical flux for SW
function pflux = f(Q)
global g
% Initialize a column vector
pflux = zeros(3,1);
h = Q(1);
if (h<=0)
    h = 0;
    u = 0;
    psi = 0;
else
    u = Q(2)/Q(1);
    psi = Q(3)/Q(1);
end

pflux(1) = h*u; % mass flux
pflux(2) = h*u^2 + 0.5 *g *h^2; % momentum flux
pflux(3) = h*u*psi;
end
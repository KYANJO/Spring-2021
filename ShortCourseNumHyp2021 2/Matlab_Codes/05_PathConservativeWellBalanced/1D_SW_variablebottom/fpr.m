function flux = fpr(Q)
global g
% Flux function
h = Q(1);
u = Q(2)/Q(1);
b = Q(3);

qs = 0;   % portata solida, =0 for fixed bottom

flux = zeros(3,1);
flux(1) = h*u;         % Linear momentum
flux(2) = h*u^2 + 1/2*g*h^2;       % Convective flux related to momentum with pressure term
flux(3) = qs;          % portata solida
end
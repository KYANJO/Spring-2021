function flux = f(Q)
% Flux function
h = Q(1);
u = Q(2)/Q(1);
b = Q(3);

qs = 0;   % solid discharge, =0 for fixed bottom

flux = zeros(3,1);

flux(1) = h*u;         % Linear momentum
flux(2) = h*u^2;       % Convective flux related to momentum (without pressure)
flux(3) = qs;          % solid discharge
end
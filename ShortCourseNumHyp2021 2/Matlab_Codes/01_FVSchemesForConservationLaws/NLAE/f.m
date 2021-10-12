function flux = f(q)
% Flux function of the conservation law

% Traffic flow equation
% u_m = 1;
% q_m = 3;
% flux = u_m*(1-q/q_m)*q;

% Burgers
flux = 0.5*q^2;

% Cubic
%  flux = 1/3*q^3;

end
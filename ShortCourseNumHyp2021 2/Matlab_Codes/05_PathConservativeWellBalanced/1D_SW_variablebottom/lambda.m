function L = lambda(Q)
% Compute eigenvalues
global g
h = Q(1);
u = Q(2)/Q(1);
c = sqrt(g*h);

L(1) = u -c;
L(2) = u;
L(3) = u+c;
end
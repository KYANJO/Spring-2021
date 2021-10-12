function BMatrix = B(Q)
global g
h = Q(1);
c2 = g*h;  % square of the celerity

BMatrix = zeros(3,3);
BMatrix(2,1) = c2;
BMatrix(2,3) = c2;
end
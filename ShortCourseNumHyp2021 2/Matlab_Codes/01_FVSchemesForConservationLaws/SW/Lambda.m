% Function to compute the eigenvalues of SW
function L = Lambda(Q)
global g
if (Q(1)<=0)
    h = 0;
    u = 0;
else
    h = Q(1);
    u = Q(2)/h;
end
c = sqrt(g*h);
L(1) = u-c;
L(2) = u;
L(3) = u+c;
end
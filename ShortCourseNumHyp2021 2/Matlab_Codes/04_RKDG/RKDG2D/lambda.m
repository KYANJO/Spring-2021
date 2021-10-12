function L= lambda(Q,n)
global U0 V0 c
if (abs(n(1))+abs(n(2))==0)
    un = sqrt(U0^2+V0^2);
else
    un = U0*n(1)+V0*n(2);
end
L = [un-c,un,un,un+c];
end
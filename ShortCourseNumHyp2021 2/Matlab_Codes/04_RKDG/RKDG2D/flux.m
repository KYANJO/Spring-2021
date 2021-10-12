function [f,g] = flux(Q)
global U0 V0 c
c2=c^2;
f =[U0*Q(1);
    U0*Q(2)+Q(4);
    U0*Q(3);
    U0*Q(4)+c2*Q(2)];
g =[V0*Q(1);
    V0*Q(2);
    V0*Q(3)+Q(4);
    V0*Q(4)+c2*Q(3)];
end
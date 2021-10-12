function f = PDEFlux(Q)
global pdetype gamma
switch(pdetype)
    case(0) % Euler
        p = (gamma-1)*(Q(3)-0.5*Q(2)^2/Q(1));
        f(1) = Q(2);
        f(2) = Q(2)^2/Q(1)+p;
        f(3) = Q(2)/Q(1)*(Q(3)+p);
end
end
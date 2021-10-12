% nonlinear flux tensor in 2D 
function flux=F(Q) 
global gamma grav eqntype
u = Q(2)/Q(1); 
v = Q(3)/Q(1);
switch(eqntype)
    case(1) % SWE 
        p = 0.5*grav*Q(1)^2;                         
        flux = [   Q(2),     Q(3); 
         u*Q(2)+p,   v*Q(2);
         u*Q(3),     v*Q(3)+p ]; 
    case(2) % Euler
        p = (gamma-1)*( Q(4) - 0.5*Q(1)*(u^2+v^2) ); 
        flux = [   Q(2),     Q(3); 
         u*Q(2)+p,   v*Q(2);
         u*Q(3),     v*Q(3)+p; 
         u*(Q(4)+p), v*(Q(4)+p) ]; 
end

end
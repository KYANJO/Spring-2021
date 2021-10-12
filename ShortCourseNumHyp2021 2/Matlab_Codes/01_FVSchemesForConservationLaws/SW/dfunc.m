% Derivata della funzione da azzerare
function y = dfunc(hs,hL,hR,uL,uR)
eps = 1e-7;
% Approsimazione della derivata utilizzando differeze finite centrale
y = (func(hs+eps,hL,hR,uL,uR) - func(hs-eps,hL,hR,uL,uR))/(2*eps);
end
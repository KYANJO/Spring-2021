% Funzione da azzerare per trovare h*
function y = func(hs,hL,hR,uL,uR)
y = phi(hs,hL) + phi(hs,hR) + uR-uL;
end
% Funzione che contine le relazioni di Rankine-Hugoniot, le invarianti di
% Riemann e le condizioni di entropia di Lax
function y = phi(hs,hLR)
global g
% Condizione di entropia di Lax
if (hs>hLR)
    % shock (condizione di Rankine-Hugoniot)
    y = sqrt(0.5*g*(hs+hLR)/(hs*hLR))*(hs-hLR);
else
    % rarefazione (invarianti di Riemann)
    y = 2*sqrt(g)*(sqrt(hs)-sqrt(hLR));
end
end
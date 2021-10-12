% Risolutore di Riemann esatto per le equazioni delle acque basse
% Input:
%   - Ql , QR : vettori di stato a sinistra e destra
%   - xi      : parametro di similitudine (xi=x/t)
% Output:
%   - Q       : soluzione del problema Q=Q(xi)

function Q = ExactRiemannSW(QL,QR,xi)
global g
hL = QL(1);                 % tirante a sinistra
hR = QR(1);                 % tirante a destra
uL = QL(2)/QL(1);           % velocità a sinistra
uR = QR(2)/QR(1);           % velocità a destra
psiL = QL(3)/QL(1);         % concentrazione a sinistra
psiR = QR(3)/QR(1);         % concentrazione a destra

hs = Newton(hL,hR,uL,uR);   % calcolare la h* con il metodo di Newton
us = uL - phi(hs,hL);       % calcolare la u* secondo la definizione

if (xi<=us)
    % a sinistra di us
    psi = psiL;
    if (hs>hL)
        % shock a sinistra
        s = uL -sqrt(0.5*g*hs/hL*(hL+hs));
        if (xi<=s)
            h = hL;
            u = uL;
        else
            h= hs;
            u = us;
        end
    else
        % rarefazione a sinistra
        head = uL - sqrt(g*hL);
        tail = us -sqrt(g*hs);
        if (xi<=head)
            % sinistra
            h = hL;
            u = uL;
        elseif (xi>=tail)
            % destra
            h = hs;
            u = us;
        else
            % dentro la rarefazione
            h = ( (uL+2*sqrt(g*hL)-xi)/3 )^2/g;
            u = xi +sqrt(g*h);
        end            
    end
else
    % a destra di us
    psi = psiR;
    if (hs>hR)
        % shock a destra
        s = uR + sqrt(0.5*g*hs/hR*(hs+hR));
        if (xi<=s)
            h = hs;
            u = us;
        else
            h = hR;
            u = uR;
        end
    else
        % rarefazione a destra
        head = uR + sqrt(g*hR);
        tail = us + sqrt(g*hs);
        if (xi>=head)
            % destra
            h = hR;
            u = uR;
        elseif (xi<=tail)
            % destra
            h = hs;
            u = us;
        else
            % dentro la rarefazione
            h = ( (xi-uR+2*sqrt(g*hR) )/3 )^2/g;
            u = xi - sqrt(g*h);
        end  
    end
end
Q = [h; h*u; h*psi]; % vettore delle variabile conservate
end
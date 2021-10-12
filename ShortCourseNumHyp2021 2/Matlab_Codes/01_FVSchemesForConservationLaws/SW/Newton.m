% Metodo di Newton per trovare h*
function hs = Newton(hL,hR,uL,uR)
global g
% Valore al primo tentativo
% hs = 0.5 *(hL+hR);    % media aritmetica
hs = ( sqrt(hL)+sqrt(hR) - (uR-uL)/2/sqrt(g) )^2/4;

tol = 1e-12;            % tolleranza
MaxNewton = 100;        % Numero massimo di iteranti permesso
for iNewton = 1:MaxNewton
    gk = func(hs,hL,hR,uL,uR);  % funzione da azzerare
    res = abs(gk);
    if (res<tol)
        break                   % soluzione trovata
    end
    dg = dfunc(hs,hL,hR,uL,uR); % derivata alla funzione
    dh = -gk / dg;              % passo di nNewton
    delta = 1;                  % fattore di scala 0<delta<1 per ridurre 
                                %   il passo di Newton
    for inner=1:20
        if (abs(func(hs+dh*delta,hL,hR,uL,uR))>= res )
            % se il residuo aumenta ridurre il passo di Newton di un
            %   fattore 2
            delta = 0.5*delta;
        else
            %se il residuo non aumenta usciamo dal ciclo
            break
        end
    end
    % aggiornamento soluzione
    hs = hs +delta*dh;
end

end
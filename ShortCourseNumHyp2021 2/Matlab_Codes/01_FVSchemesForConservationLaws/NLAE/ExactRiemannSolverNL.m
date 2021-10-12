function q = ExactRiemannSolverNL(qL,qR,xi)
% Funtion for the computation of the exact solution of the Riemann Problem
% for q_t + f_x = 0
% Lax entropy condition
if (a(qL)>a(qR))
    % Shock wave
    s = (f(qR)-f(qL)) /(qR-qL);
    if (xi > s) % at the right of the shock wave
        q = qR;
    else
        q = qL; % at the left of the shock wave
    end
else
    % Rarefaction wave
    if (xi <= a(qL)) % at the left of the rarefaction wave
        q = qL;
    elseif (xi >= a(qR)) % at the right of the rarefaction fan
        q = qR;
    else
        q = Newton(xi); % Compute the solution of g(q) = a(q) -xi = 0
    end  
end

end
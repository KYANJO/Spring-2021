function q = Newton(xi)
% Funtion to solve the scalar non linear equation g(q)=0 using Newton
% method
% Initial guess
q = xi;
% Tolerance
tol = 1.e-12;
% Maximum number of iterants
MaxNewton = 100;
% Newton loop
for i = 1:MaxNewton
    res = abs(g(q,xi));
    if (res < tol)
        break
    end    
    % Compute Newton step
    dq = -g(q,xi)/dg(q,xi);
    % Update solution
    q = q+dq;
end

end
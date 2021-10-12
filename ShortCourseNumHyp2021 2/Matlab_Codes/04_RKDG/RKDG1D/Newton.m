function q=Newton(xi) 
tol = 1e-12;        % tolerance 
MaxNewton = 100;    % max. number of Newton iterations 
q = xi;             % initial guess 
for iNewton=1:MaxNewton
    gk = g(q,xi);   % compute the nonlinear function 
    if(abs(gk)<tol)
        break
    end
    dgk = dg(q,xi); % compute the derivative 
    dq  = -gk/dgk;  % compute the Newton step 
    q   = q + dq;   % update solution 
end

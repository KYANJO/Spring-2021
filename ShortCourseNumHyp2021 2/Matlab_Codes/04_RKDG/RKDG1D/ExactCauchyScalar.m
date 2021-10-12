% exact solution of the Cauchy problem of the nonlinear scalar conservation law 
%   q_t + f_x = 0   (PDE) 
%   q(x,0) = h(x)   (IC) 
% making use of the method of characteristics 
% Be careful: this code can NOT solve problems with shocks, but only smooth problems
% where characteristics do not intersect or have not intersected yet !! 
function q=ExactCauchyScalar(x,t) 
x0 = x-t*a(h(x));   % initial guess
tol = 1e-12;        % tolerance 
MaxNewton = 100;    % max. number of Newton iterations 
% Newton iterations 
for iNewton=1:MaxNewton
    res = abs(gs(x0,x,t)); 
    if(res<tol)
        break
    end
    dx0 = -gs(x0,x,t)/dgs(x0,x,t);  % Newton step 
    x0  = x0 + dx0;               % update the solution 
end
% The exact solution is given by the initial condition at the foot of the characteristics 
q = h(x0); 
    
    
% Exact solution of the Riemann problem for scalar conservation laws 
function q=ExactRiemannScalar(qL,qR,xi)
% Lax entropy condition 
if(a(qL)>a(qR)) 
    % shock
    s = (f(qR)-f(qL))/(qR-qL); 
    if(xi<=s)
        q = qL;
    else
        q = qR;
    end    
else
    % rarefaction 
    if(xi<=a(qL))
        q = qL;
    elseif(xi>=a(qR))
        q = qR;
    else
        q = Newton(xi); 
    end    
end
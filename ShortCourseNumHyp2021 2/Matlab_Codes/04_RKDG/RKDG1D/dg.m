function dgdq=dg(q,xi)
% central finite difference approximation for dg/dq 
eps = 1e-7; 
dgdq = ( g(q+eps,xi)-g(q-eps,xi) )/(2*eps); 
% approximation of the derivative of the function gs using a 
% central finite difference 
function y=dgs(x0,x,t) 
eps = 1e-7; 
y = ( gs(x0+eps,x,t)-gs(x0-eps,x,t) )/(2*eps); 


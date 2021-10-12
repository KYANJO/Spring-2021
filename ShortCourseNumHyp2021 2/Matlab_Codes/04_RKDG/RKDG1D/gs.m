% function whose root we have to find in order to get the exact solution of
% the Cauchy problem of the scalar conservarion law 
function y=gs(x0,x,t) 
y = x - x0 - a(h(x0))*t;

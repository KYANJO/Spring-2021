function y = dg(q,xi)
% Approximate derivative
eps = 1.e-7;
y = (g(q+eps,xi)-g(q-eps,xi))/(2*eps);
end
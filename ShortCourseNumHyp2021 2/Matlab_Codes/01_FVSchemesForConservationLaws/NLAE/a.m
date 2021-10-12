function dfdq = a(q)
% Characteristic velocity a(q)=df/dq
% for the scalar non linear conservation law
% q_t + f_x = 0

% Traffic flow equation
% u_m = 1;
% q_m = 3;
% dfdq = u_m*(1-2*q/q_m);

% Burgers
 dfdq = q;

% Cubic
% dfdq = q^2;
end
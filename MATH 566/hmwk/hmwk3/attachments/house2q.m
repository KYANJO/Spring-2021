%==========================================================================
% function computes tcomputes the corresponding m-by-m orthogonal matrix Q.
% input : real m by n matrix V
% output : real m by m matrix Q
%==============================================

function Q = house2q(V)

[m,n] = size(V);  

x = eye(m);

for k = n:-1:1
    x(k:m,:) = x(k:m,:) - 2*V(k:m,k)*(V(k:m,k)'*x(k:m,:));
end
Q = x;
end
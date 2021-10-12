clear all;
close all;
%m equally spaced points over [0,1]
m = 50; n=12;

% Vandermonde matrix t
t = zeros(m,n);
for i = 1:n
    for j = 1:m
        t(j,i) = ((j-1)/(m-1))^(n-i);
    end
end

%fliping the vandermonde matrix t to form A
A = fliplr(t);

%fuction f
tj = zeros(m,1);
for j = 1:m
    tj(j) = (j-1)/(m-1);
end

f = cos(4*tj);

format long
%block matrix Ablock
Ablock = [eye(m) A; A' zeros(n)];
b = [f;zeros(n,1)];
[qb,rb] = qr(Ablock); xb = rb\(qb'*b);

%solution x, extracted from xb
xblock = xb(m+1:end)

%solution r
fprintf('Solution of r\n\n');
r = xb(1:m);

%(a). normal equations
x = (A'*A)\(A'*f);


%(e). QR decomposition using inbuilt Householder
[q,r] = qr(A); xh = r\(q'*f);

Table = table(x,xh,xblock, 'VariableNames',{'Normal equation','Builtin function','Block Matrix'})

residual = norm((xh - xblock),2);

fprintf('The solution, xblock, solution is almost similar to the computed solutions although some little values are off. \n And the residual %f is small. \n',residual);

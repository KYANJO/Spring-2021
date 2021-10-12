close all;
clear all;

%Complete Pivoting
fprintf('no4b.\n\n')
n = 12;
A = vandermonde(n);

b = A(:,n);

[LU,p,q,gf,L,U]= lucp(A);
P = eye(n); P=P(p,:);
Q = eye(n); Q = Q(:,q);
yc = forsub(L,P*b);
zc = backsub(U,yc);

fprintf('Compute pivoting solution \n')
xc = Q*zc

x = A\b
residual = norm(b - A*xc);
fprintf('The value of the max-norm of the residual = %f\n\n',residual);

fprintf('no4c.\n\n')
%Gaussian Elimination with Partial pivoting
[Lp,Up,P1,g] = lupp(A);
yp = forsub(Lp,P1*b);

fprintf('Partial pivoting solution \n')
xp = backsub(Up,yp)

residual_xp = norm(b - A*xp);
fprintf('The value of the max-norm of the residual = %f\n\n',residual_xp);

fprintf('Results from part (b), are exactly the same as the exact solution x, since even the residual\n is zero. But for part(c), the solution doesnot properly approximate the actual solution, according to the residual computed. \n\n')

function A = vandermonde(n)
    
    t = zeros(n,n);
    for i = 1:n
        for j = 1:n
            t(j,i) = j^(n-i);
        end
    end

    %fliping the vandermonde matrix t to form A
    A = fliplr(t);

end


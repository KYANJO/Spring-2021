% Generalize Minimum Residual (GMRES) method for solving Ax = b 
clear 
close all

%dimensions of A 
m = 4;
%matrix A
A = 2*eye(m) + 0.5*randn(m)/sqrt(m); 
%vector b
b = ones(m,1);

X = GMRES(A,b)

x = gmres(A,b)

function X = GMRES(A,b)
    m = length(A);
    q = zeros(m,m+1);
    %x = zeros(m,1);
    %r = zeros(m,1);
    h = zeros(m+1,m);
    %H = zeros(m+1,m);
    e1 = zeros(m+1,1); e1(1) = 1;
    em = zeros(m,1); em(m) = 1;
    nb = norm(b,2);
    q(:,1)= b/nb;
    tol = 1e-10;
    for n = 1:m
        v = A*q(:,n);
        for j = 1:n
           h(j,n) = q(:,j)'*v;
           v = v - h(j,n)*q(:,j);
        end
        h(n+1,n) = norm(v,2);
        q(:,n+1) = v/h(n+1,n);
        H = triu(h,-1);
        %H = [H; h(n+1,n)*em']
        %H(1:m,1:m) = Hn; 
        %H(m+1,m)= h(n+1,n) 
        b1=nb*e1;
        %b1 =[b1;0] 
        y = H\b1;
        x = q.*y;
        %xn = x(n);
        %calculated the residual
        r = A*x - b;
  
        if norm(r,2)/nb <= tol
            break
        end
    end
    X = x;
end
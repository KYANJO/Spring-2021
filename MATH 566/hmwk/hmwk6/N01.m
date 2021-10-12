% Generalize Minimum Residual (GMRES) method for solving Ax = b 
clear 
close all

%dimensions of A 
m = 200;
%matrix A
A = 2*eye(m) + 0.5*randn(m)/sqrt(m); 
%vector b
b = ones(m,1);
%tolerance
tol = 1e-10;
X = GMRES(A,b,tol);

%x = gmres(A,b);

function X = GMRES(A,b,tol)
    m = length(A);
    q = zeros(m,m);
    h = zeros(m,m);
    nb = norm(b,2);
    q(:,1)= b/nb;
    
    for n = 1:m
        v = A*q(:,n);
        for j = 1:n
           h(j,n) = q(:,j)'*v;
           v = v - h(j,n)*q(:,j);
        end
        h(n+1,n) = norm(v,2);
        q(:,n+1) = v/h(n+1,n);
        H = h(1:n+1,1:n);
        
        b1=nb*speye(n+1,1); 
        y = H\b1;
        xn = q(:,1:n)*y;
        
        %calculated the residual
        r = A*xn - b;
        
       %relative residual
        R = norm(r,2)/nb;
        
        iter = n;      
        if (norm(r,2) <nb * tol)
            break
        end
          
    end
    fprintf("The code converged at %d iterations to solution with relative residual %e\n",iter,R);

    X = xn;
end

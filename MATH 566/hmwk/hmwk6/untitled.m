% Generalize Minimum Residual (GMRES) method for solving Ax = b 
clear 
close all

%dimensions of A 
m = 200;
%matrix A
A = 2*eye(m) + 0.5*randn(m)/sqrt(m); 
%vector b
b = ones(m,1);


    m = length(A);
    q = zeros(m,m+1);
    h = zeros(m+1,m);
    x = zeros(m,1);
    r = zeros(m,1);
    e1 = zeros(m+1,1); e1(1) = 1;
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
        
        H = h(1:n+1,1:n);
        
         b1=nb*e1(1:n+1);
         y = H\b1;
         x(1:n) = q(1:n,1:n)*y;
        
%          %calculated the residual
          r = A*x - b;
%    
%          if norm(r,2)/nb <= tol
%              break
%          end
%         %relative residual
%          R = norm(r,2)/nb;
%          iter = n;        
      end
%     fprintf("The code converged at %d iterations to solution with relative residual %e\n",iter,R);









%X = GMRES(A,b);

%x = gmres(A,b);

function X = GMRES(A,b)
    m = length(A);
    q = zeros(m,m+1);
    h = zeros(m+1,m);
    e1 = zeros(m,1); e1(1) = 1;
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
        %H = triu(h,-1)
        
        b1=nb*e1(1:n); 
        y = H\b1;
        x = q(1:m,1:m)*y
        
        %calculated the residual
        r = A*x - b;
  
        if norm(r,2)/nb <= tol
            break
        end
        %relative residual
        R = norm(r,2)/nb;
        iter = n;        
    end
    fprintf("The code converged at %d iterations to solution with relative residual %e\n",iter,R);

    X = x;
end

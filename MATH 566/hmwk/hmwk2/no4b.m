% Program that implements the modified Gram-Schmidt algorithm for
% computing the QR factorization of a m by n matrix A, where m>=n.

function [q,r] = MGS(A)

    [m,n] = size(A);
    
    v = zeros(m,n);
    r = zeros(n,n);
    q = zeros(m,n);
    
    for i = 1:n
        v(:,i) = A(:,i);
    end
    
    for i = 1:n
       r(i,i) = norm(v(:,i),2);
       q(:,i) = v(:,i)/r(i,i);
       for j = i+1:n
          r(i,j) = q(:,i)'*v(:,j);
          v(:,j) = v(:,j) - r(i,j)*q(:,i);
       end
       
    end

end
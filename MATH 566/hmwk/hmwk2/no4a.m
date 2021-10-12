% Program that implements the classical Gram-Schmidt algorithm for
% computing the QR factorization of a m by n matrix A, where m>=n.

%m equally spaced points over [0,1]

%Classical Gram-Schmidt 
function [q,r] = CGS(A)

[m,n] = size(A);

q = zeros(m,n);
r = zeros(n,n);
v = zeros(m,n);

for j=1:n
    
   v(:,j) = A(:,j);
   
   for i = 1:j-1
       
      r(i,j) = q(:,i)'*A(:,j);
      v(:,j) = v(:,j) - r(i,j)*q(:,i); 
      
   end
   
   r(j,j) = norm(v(:,j),2);
   q(:,j) = v(:,j)/r(j,j);
   
end
end
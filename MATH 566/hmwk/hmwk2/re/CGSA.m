% ==============================================================
% This function implement the Classical Gram-Schmidt algorithm
%
% Input: 
%       Matrix A
%
% Outputs: 
%       Matrix Q
%       Matrix R, an upper triangular matrix
% ==============================================================

function [Q, R] = CGSA(A)

[m,n] = size(A);

Q = zeros(m,n);
R = zeros(n,n);
v = zeros(m,n);

for j = 1:n
    
    v(:,j) = A(:,j);
    
    for i = 1:j-1
        
        R(i,j) = Q(:,i)'*A(:,j);
        v(:,j) = v(:,j) - R(i,j)*Q(:,i);
    
    end
    
    R(j,j) = norm(v(:,j),2);
    Q(:,j) = v(:,j)/R(j,j);
   
    
end
end
   
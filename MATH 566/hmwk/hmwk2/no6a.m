

function [V,R] = house(A)

[m,n] = size(A);

R = A;

V = zeros(m,n);

for k=1:n
    
   x = R(k:m,k); 
   V(k:m,k) = sign(x(1))*(norm(x,2))*eye(m-k+1,1) + x;
   V(k:m,k) = V(k:m,k)/(norm(V(k:m,k),2));
   R(k:m,k:n) = R(k:m,k:n) - 2*V(k:m,k)*(V(k:m,k)'*R(k:m,k:n));
   
end

end
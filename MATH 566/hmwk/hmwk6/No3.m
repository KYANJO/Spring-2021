fprintf("Compute an eigenvalue–eigenvector pair of the matrix\n\n");
clear
close all;
%initial vector
x0 = [1 1 1 1 1 1]'./(6^0.5);
%matrix A
m = 6; n = 6;
A = matrix(m,n)

%tolerance
ep = 1e-10;

%fuction rqi
[v,lam] = rqi(A,x0,ep)
fprintf("Eigen value and its corresponding eigen vector to which the code converges are v and lam.\n\n");

%verification
[V,D] = eig(A);

fprintf("The error in the approximate eigen value.\n\n");
error_lam = abs(D(4,4) - lam)
fprintf("The error in the approximate eigen value.\n\n");
error_v = abs(V(:,4) - (-v))
fprintf("Hence v and lam are well aprroximated since the error is too small.\n\n");

function A = matrix(m,n)
    A = zeros(m,n); 
    for i = 1:n
       for j = 1:m
          if i ==j 
             A(i,j) = -2;
          elseif i == j+1
             A(i,j) = 1;
          elseif i == j-1
             A(i,j) = 1;
          end
          A(1,2) = 2;
       end

    end
end
clear all
close all

fprintf('4b). Use your code to compute the growth factor of the 21 by 21 matrix. \n\n')

m = 21;
A = zeros(m,m);
for i = 1:m
   for j = 1:m
       if i == j
            A(i,j) = 1;
       elseif i > j
           A(i,j) = -1;
       else
           A(i,m) = 1;
       end
   end
end

[L,U,P,g] = lupp(A);

fprintf('The growth factor for A is: %d\n',g);



clear all

close all

clc
 
 
a1 = [1 2 3 4 5 6 7 8 7];
a2 = [2 4 6 5 2 8 3 5 1];

a3 = [2 3 1 2 6 9 4 5 6];

A = diag(a1) + fliplr(diag(a2));

a4 = [1 3 4 6 7 9 2 7 8];
a5 = [4 1 2 3 9 7 5 4 2];
a6 = [6 7 3 4 9 1 9 4 5];

A(:,1) = a3;
A(1,:) = a4;
A(9,:) = a6;
A(:,9) = a5;
  
x_exact = [1 2 5 8 9 12 18 5 4];

b = A*x_exact';

[a,x] = Gauss(A,b);

 function [A,x] = Gauss(A,b)

     A = [A,b];

     [m,n] = size(A);



    for i = 2:m-1
      if(m-i+1~=i)
          A(m-i+1,:) = A(m-i+1,:) - A(i,:)*A(m-i+1,i)/A(i,i);
          
      end

   end

    for k = 1:m-1

        A(m,:) = A(m,:) - A(k,:)*A(m,k)/A(k,k);
        
    end

    A1 = A(1,:);
    A(1,:) = A(m,:);

    A(m,:) = A1;


    for j = 2:m

        A(j,:) = A(j,:) - A(1,:)*A(j,1)/A(1,1);
    end

    for k = 2:m-1

        A(m,:) = A(m,:) - A(k,:)*A(m,k)/A(k,k);
        
    end
    
    
    x = zeros(1,m);

    

    for s = m:-1:1

       c = 0;

       for k = s:m

           c = c + A(s,k)*x(k);

       end

       x(s) = (A(s,n)-c)/A(s,s);

    end
     
 end
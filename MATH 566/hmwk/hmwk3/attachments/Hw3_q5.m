clear all
close all
clc
% Dimension on the matrix
   

% Vendermond matrix

cdnb = zeros(29,1);

for n = 2:30
    
    m = n;
    
    A = vendermond(m,n);
    
    cd =cond(A);
    cdnb(n-1) = cd;
    
end    

cdnb1 = zeros(29,1);
for n = 2:30
    
    m = 2*n-1;
    
    A = vendermond(m,n);
    
    cd =cond(A);
    cdnb1(n-1) = cd;
    
end    

N = 2:30;
semilogy(N,cdnb,'r','linewidth',2)
%set(gca,'FontSize',13)
hold on
semilogy(N,cdnb1,'b--','linewidth',2)
xlabel('n')
ylabel('cond(A)')
title('Plot of function condition number of A')
legend('m = n', 'm = 2n-1','Location','northwest')
set(gca,'FontSize',13)


function A = vendermond(m,n)

    A = zeros(m,n);
    J = 1:m;
    t = (J-1)/(m-1);

    for j = 1:n

        A(:,j) = t.^(j-1);
    end
    
end


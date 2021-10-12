clear all
close all
clc
% Dimension on the matrix
m = 50;   
n = 12;

% Vendermond matrix
A = zeros(m,n);
I = eye(n);
J = 1:m;
t = (J-1)/(m-1);

for j = 1:n
    
    A(:,j) = t.^(j-1);
end

f = @(t) cos(4.*t);

b = f(t)';

%fprintf("(a) Use of the normal equations\n");
x = (A'*A)\(A'*b);
%disp(x)

%fprintf("(b) QR classical Gram-Schmidt\n");
[Q1, R1] = CGS(A);
y1 = Q1'*b;
%x1 = backward(R1,y1);  
x1=R1\y1;
%disp(x1)

%fprintf("(c) QR modified Gram-Schmidt\n");
[Q2, R2] = MGS(A);
y2 = Q2'*b;
%x2 = backward(R2,y2);
x2=R2\y2;
%disp(x2)

%fprintf("(d) QR Householder matrice\n");
[V,R3] = house(A);
Q3 = house2q(V);
y3 = Q3'*b;
%x3 = backward(R3,y3);%
x3 = R3\y3;
%disp(x3)

%fprintf("(e) QR builtin function\n");
[Q4,R4] = qr(A);
y4 = Q4'*b;
%x4 = backward(R4,y4);%
x4 = R4\y4;
%disp(x4)

%fprintf("(f) SVD builtin function\n");
[U,S,V]=svd(A,0);
x5= V*((U'*b)./diag(S));
%disp(x5)

Tab = table(round(x,15),round(x1,15),round(x2,15),round(x3,15),round(x4,15),...
round(x5,15),'VariableNames',{'\ Matlab','CGS','MGS','QR Householder',...
'QR Matlab','SVD'});

disp(Tab)



pa = zeros(m,1);
pe = zeros(m,1);

K = 1:n;

for i = 1:m
    
    ti = t(i).^(K-1);
    pa(i) = sum(x.*ti');
    pe(i) = sum(x4.*ti');
    
end
        

difpaf = b - pa;
difpef = b - pe;

plot(t, difpaf, 'linewidth',2)
hold on
plot(t, difpef, 'linewidth',2)
xlabel('t')
ylabel('f(t)-p(t)')
title('Difference  between the approximat polynomial in (a) and (e) and f(t)')
legend('f(t)-p(t) in (a)', 'f(t)-p(t) in (e)', 'Location', 'southeast')
set(gca,'FontSize',13)


        
        
        
        
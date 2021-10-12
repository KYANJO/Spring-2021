clear all;
close all;
%m equally spaced points over [0,1]
m = 50; n=12;

% Vandermonde matrix t
t = zeros(m,n);
for i = 1:n
    for j = 1:m
        t(j,i) = ((j-1)/(m-1))^(n-i);
    end
end

%fliping the vandermonde matrix t to form A
A = fliplr(t);

%fuction f
tj = zeros(m,1);
for j = 1:m
    tj(j) = (j-1)/(m-1);
end

f = cos(4*tj);

format long
%(a). normal equations
x = (A'*A)\(A'*f);

%(b). QR decomposition using CGS
[q_c,r_c] = CGS(A); xc = r_c\(q_c'*f);

%(c). QR decomposition using MGS
[q_m,r_m] = MGS(A); xm = r_m\(q_m'*f);

%(d). QR decomposition using Householder
[v_h,r_h] = house(A); q_h = house2q(v_h);
x_h = r_h\(q_h'*f);

%(e). QR decomposition using inbuilt Householder
[q,r] = qr(A); xh = r\(q'*f);

%(f). QR decomposition using inbuilt svd
[u,s,v] = svd(A); xs = (u*s*v')\f;

Table = table(x,xc,xm,x_h,xh,xs, 'VariableNames',{'Normal equation','CGS','MGS','Householder','Builtin function','SVD'})

%Differences and Similarities
fprintf('The Normal equation and the MGS, slightly give the same results different from SVD, CGS, built in function and the Householder, however the SVD,\n built in function and the Householder matrix give almost similar results different from CSG.\n')
%Plot the difference between AX - b
%a) the Equations method
plot(tj,(f - A*x),'-*')
hold on
%e) the Inbulit in Householder
plot(tj,(f - A*xh),'-o')
title('Ax - f against t')
xlabel('tj');ylabel('Ax - f')
legend('Equations method','Householder')



%matrix A
A = [1 2 3 ; 4 5 6 ; 7 8 7 ; 4 2 3 ; 4 2 2];
 
%calling function house
[V,R] = house(A)

%calling function house2q
Q = house2q(V)

%function qr from matlab
[q ,r] = qr(A)

fprintf('The Q and R, produced by the code satisy A = QR \n');
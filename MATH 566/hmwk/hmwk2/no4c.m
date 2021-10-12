clear all;
%m equally spaced points over [0,1]
m = 100; n=15;

% Vandermonde matrix B
t = zeros(m,n);
for i = 1:n
    for j = 1:m
        t(j,i) = ((j-1)/(m-1))^(n-i);
    end
end

%fliping the vandermonde matrix t to form A
A = fliplr(t);

%calling the function for CGS
[q_c,r_c] = CGS(A);

%infinity norm of A - QR for CGS
Nc = norm((A - q_c*r_c),inf)
Nc1 = norm((eye(n) - q_c'*q_c),inf)

%calling the function for MGS
[q_m,r_m] = MGS(A);

%infinity norm of A - QR for MGS
Nm = norm((A - q_m*r_m),inf)
Nm2 = norm((eye(n) - q_m'*q_m),inf)

fprintf("The infinity norms (A - QR) for the two methods are very small to almost zero, \n while for (I - Q'Q), for the modified algorithm its near to zero, but its \n strange for the classical one, because its 4.6752 and this is a very big \n value sice we are almost subtracting the same things. \n");

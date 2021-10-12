% Compute the weighted least squares solution using the diagonal Gaussian
% weight with t = 1/23
clear all
close all

tex = 1/23;
delta = 1;

%exact approximation
fex = cos(4*tex);

%Guassian function
w = @(t,tj,delta) exp(-(abs(t - tj)/delta).^2);

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

%Compute the weighted least square Using the Diagonal Gaussian  weight, W
W = diag(w(tex,tj,delta));

format long
%Report the polynomial coefficients of the weighted least squares solution.
fprintf('Polynomial coefficients of the weighted least squares solution \n');
[qw,rw] = qr(W*A); xw = rw\(qw'*(W*f))

%Non Weighted least squares solution xh
fprintf('Polynomial coefficeients non weighted least squares solution \n');
[q,r] = qr(A); xh = r\(q'*f)

%Report the value of the polynomial with these coefficients at t =1/23
% Vandermonde matrix t
tc = zeros(m,n);
for i = 1:n
    for j = 1:m
        tc(j,i) = (1/23)^(n-i);
    end
end

%fliping the vandermonde matrix t to form A
Ac = fliplr(tc); 

%value of the polynomial at t =1/23
pw = Ac*xw; pw(11);
fprintf('Polynomial value computed using weighted coeffieients:')
disp(pw(11));

%Compare these coefficients 
%for non weighted coefficients
pnonw = Ac*xh; pnonw(11);
fprintf('Polynomial value computed using non-weighted coeffieients:');
disp(pnonw(11));
fprintf('Exact Polynomial value computed directly:')
disp(fex)
%Which method provides better approximation?
fprintf('Comparing the three polynomial values, its clear that the onw computed with the weighted \n coefficients best approximates the polynomial compared to the one computed with non-weighted coefficients.\n')

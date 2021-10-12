function x = dct(x,dim)
%  DCT   Discrete Cosine Transform computed using the fast Fourier Transform.
%     X = dct(x) computes the Discrete Cosine Transform (DCT) of the columns of X.
%
%     X_k = 1/N*(1/2*x_0 + \sum_{j=1}^{N-1} x_j*cos(\pi*j*k/N) + 1/2*(-1)^k*x_N)
%
%     k = 0,...,N, where N is the column length of X.
%
%     X = dct(x,dim) computes the DCT along the dimension specified.
%     if dim = 1 (default) then the DCT is along the columns.
%     if dim = 2 then the DCT is along the rows.
%
%  See also idct.

[m,n] = size(x);

if nargin == 1
    dim = 1;
end

if dim == 1
    xe = real( fft( [x; x(m-1:-1:2,:)]/(2*(m-1)),[],dim ) );
    x = xe(1:m,:);
elseif dim == 2
    xe = real( fft( [x x(:,n-1:-1:2)]/(2*(n-1)),[],dim ) );
    x = xe(:,1:n);
else
    error('dct:dimUnknown','DCT dimension not available, select 1 or 2');
end

end


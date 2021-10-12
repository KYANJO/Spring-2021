function x = idct(x,dim)
%  IDCT   Inverse Discrete Cosine Transform computed using the fast Fourier Transform.
%     x = idct(X) computes the inverse Discrete Cosine Transform (DCT) of the columns of X.
%
%     x_j = X_0 + 2\sum_{k=1}^{N-1} X_k*cos(\pi*j*k/N) + 1/2*(-1)^j*X_N,
%
%     j = 0,...,N, where N is the column length of X.
%
%     X = idct(x,dim) computes the inverse DCT along the dimension specified.
%     if dim = 1 (default) then the inverse DCT is along the columns.
%     if dim = 2 then the inverse DCT is along the rows.
%
%  See also dct.
if nargin == 1
    dim = 1;
end

[m,n] = size(x);

if dim == 1
    scale = (2*m-2);
elseif dim == 2
    scale = (2*n-2);
else
    error('idct:dimUnknown','IDCT dimension not available, select 1 or 2');
end
x = scale*dct(x,dim);

end


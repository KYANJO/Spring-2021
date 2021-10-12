function x = dst(x,dim)
%  DST   Discrete Sine Transform computed using the fast Fourier Transform.
%     X = dst(x) computes the Discrete Sine Transform (DST) of the columns of X.
%
%     X_k = 1/N*\sum_{j=1}^{N-1} x_j*sin(\pi*j*k/N), k=1,...,N-1,
%
%     where N-1 is the column length of x.  Leading zeros should not
%     be present input
%
%     X = dst(x,dim) computes the DST along the dimension specified.
%     if dim = 1 (default) then the DST is along the columns.
%     if dim = 2 then the DST is along the rows.
%  
%     See also idst.

if nargin == 1
    dim = 1;
end

[m,n] = size(x);

if dim == 1
    xe = zeros(2*m+2,n);
    xe(2:m+1,:) = x;
    xe = imag(fft(xe));
    x = 1/(m+1)*xe(2:m+1,:);
elseif dim == 2
    xe = zeros(m,2*n+2);
    xe(:,2:n+1) = x;
    xe = imag(fft(xe,[],2));
    x = 1/(n+1)*xe(:,2:n+1);
else
    error('dst:dimUnknown','DST dimension not available, select 1 or 2');
end

end


function x = idst(x,dim)
%  IDST   Inverse Discrete Sine Transform computed using the fast Fourier Transform.
%     x = idst(X) computes the inverse Discrete Sine Transform (DST) of the columns of X.
%
%     x_j = 2\sum_{k=1}^{N-1} X_k*sin(\pi*j*k/N), j=1,...,N-1,
%
%     where N-1 is the column length of X.  Leading zeros should not
%     be present input.
%  
%     X = idst(x,dim) computes the inverse DST along the dimension specified.
%     if dim = 1 (default) then the inverse DST is along the columns.
%     if dim = 2 then the inverse DST is along the rows.
%
%     See also dst.
if nargin == 1
    dim = 1;
end

[m,n] = size(x);

if dim == 1
    scale = 2*(m+1);
elseif dim == 2
    scale = 2*(n+1);
else
    error('idst:dimUnknown','IDST dimension not available, select 1 or 2');
end
x = scale*dst(x,dim);


end


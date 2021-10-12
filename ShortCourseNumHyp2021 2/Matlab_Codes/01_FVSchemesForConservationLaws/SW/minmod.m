% Function for the computation of the slope limiter (minmod)
function c = minmod(a,b)
% Compute input variables length
nVar = length(a);
% Allocate memory for a column vector
c = zeros(nVar,1);
for i = 1:nVar
    if (a(i)*b(i)<= 0 )
        c(i) = 0;
    else
        if (abs(a(i))<abs(b(i)))
            c(i) = a(i);
        else
            c(i) = b(i);
        end
    end
end
end
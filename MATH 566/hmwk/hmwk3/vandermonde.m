

function A = vandermonde(m,n)
    
    t = zeros(m,n);
    for i = 1:n
        for j = 1:m
            t(j,i) = ((j-1)/(m-1))^(n-i);
        end
    end

    %fliping the vandermonde matrix t to form A
    A = fliplr(t);

end
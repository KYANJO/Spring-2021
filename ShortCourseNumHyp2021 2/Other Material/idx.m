function j=idx(i)
global IMAX
j=i; 
if(j<=0)
    j=j+IMAX;
elseif(j>=IMAX+1)
    j=j-IMAX; 
end


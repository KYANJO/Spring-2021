function j = idx(i)
% Identify the index of neighbour boundaries when considering PBC
global IMAX
j = i;
if (j<=0)
    j = j+IMAX;
elseif (j>=IMAX+1)
    j = j-IMAX;
end
end
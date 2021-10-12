% The minmod limiter function. 
%  Syntax:
%    limslope=minmod(slope1,slope2)
function limslope=minmod(slope1,slope2)
  omega    = abs(slope1)<abs(slope2);
  limslope = omega.*slope1 + (1-omega).*slope2;
  signch   = slope1.*slope2>0;
  limslope = limslope.*signch;
end
function w_l = interpolate(w_lm1, parents)
% Interpolation from coarse to fine.
% Uses the algorithm from Jung/Langer: Methode der finiten Elemente
% fuer Ingenieure, Teubner, 2001, p. 279.
%
% INPUT:   w_lm1    vector on coarse grid
%          parents  father-son relation
%
% OUTPUT:  w_l      vector on fine grid
%
% VERSION 1.0
% DATE 25.3.2004
% EMAIL bernd@flemisch.net
  
  n_lm1 = length(w_lm1);
  n_l = n_lm1 + size(parents, 2);
  w_l = zeros(n_l,1);
  w_l(1:n_lm1) = w_lm1;
  for i = n_lm1+1:n_l
    idx1 = parents(1, i-n_lm1);
    idx2 = parents(2, i-n_lm1);
    w_l(i) = 0.5*(w_l(idx1) + w_l(idx2));
  end
  
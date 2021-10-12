function d_lm1 = restrict(d_l, parents, isdir)
% Restriction from fine to coarse.
% Uses the algorithm from Jung/Langer: Methode der finiten Elemente
% fuer Ingenieure, Teubner, 2001, p. 279.
%
% INPUT:   d_l      vector on fine grid
%          parents  father-son relation
%
% OUTPUT:  d_lm1    vector on coarse grid
%
% VERSION 1.0
% DATE 25.3.2004
% EMAIL bernd@flemisch.net

  n_l = length(d_l);
  n_lm1 = n_l - size(parents, 2);
  d_lm1 = d_l(1:n_lm1);
  for i = n_lm1+1:n_l
    idx1 = parents(1, i-n_lm1);
    idx2 = parents(2, i-n_lm1);
    d_lm1(idx1) = d_lm1(idx1) + 0.5*d_l(i);
    d_lm1(idx2) = d_lm1(idx2) + 0.5*d_l(i);
  end
  [i, j] = find(isdir(1:n_lm1) == 1);
  d_lm1(j) = 0;
  
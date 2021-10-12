function edges = get_edges(p, t)
% Locates the edges of a finite element mesh.
% Especially useful for the implementation of transfer operators in
% the context of multigrid methods, since the endpoints of an edge are
% exactly the two fathers of a newly created FE node.
%
% INPUT:   p  2 by np matrix containing the coordinates of the
%             finite element nodes
%          t  4 by nt matrix containing the indices (with
%             respect to p) of the vertices of the triangulation
%
% OUTPUT:  edges  2 by ne matrix, ne the number of edges. The i-th
%                 column contains the indices (with respect to p) of the
%                 edge endpoints
%
% VERSION 1.0
% DATE 25.3.2004
% EMAIL bernd@flemisch.net

  np = size(p,2);
  nt = size(t,2);
  it = (1:nt); 
  ip1 = t(1,it);
  ip2 = t(2,it);
  ip3 = t(3,it);
  A = sparse(ip1, ip2, -1, np, np);
  A = A + sparse(ip2, ip3, -1, np, np);
  A = A + sparse(ip3, ip1, -1, np, np);
  A = -((A + A.') < 0);
  [i1,i2] = find(A == -1 & A.' == -1);
  i = find(i2 > i1);
  i1 = i1(i)';
  i2 = i2(i)';
  edges = [i1; i2];

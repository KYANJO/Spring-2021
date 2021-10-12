function isdir = dirichlet_flags(p, e, b)
% set dirichlet flags for the nodes of the finite element mesh
%
%  INPUT:   p   matrix containing the coordinates of the FE nodes
%           e   matrix containing the information of the boundary edges
%           b   boundary condition data, either a string with the
%               name of a boundary M-file (see Matlab help to
%               pdebound for details) or a boundary condition
%               matrix created by the PDE toolbox (see Matlab help
%               to assemb for details)
%               
%
%  OUTPUT:  isdir  vector of length np, the i-th entry corresponds
%                  to the i-th FE node, 1 for a Dirichlet boundary
%                  node, 0 for an inner node or a Neumann boundary
%                  node 
%
%  VERSION 1.0
%  DATE 25.3.2004
%  EMAIL bernd@flemisch.net

  isdir = zeros(1, size(p,2));
  ne = size(e, 2);
  if (ischar(b))
    [q,g,h,r] = feval(b,p,e,[],[]);
    for k = 1:ne
      if (h(k) ~= 0)
	isdir(e(1,k)) = 1;
      end
    end
    for k = ne+1:2*ne
      if (h(k) ~= 0)
	isdir(e(2,k-ne)) = 1;
      end
    end
  else
    for k = 1:ne
      if (b(2, e(5,k)) == 1)
	isdir(e(1:2,k)) = [1 1];
      end
    end
  end
  
      
      
      
      
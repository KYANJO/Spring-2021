% Rayleigh quotient iteration
%==========================================================================
% function rqi is a powerful method for finding an eigenvalue–eigenvector 
% pairs of certain matrices (especially symmetric tridiagonal ones!
% input: A - matrix
%        x0 - intial starting
%        ep - tolerance
% output: v and lam are the eigenvalue-eigenvector pair that the algorithm
%         converged to.
%==========================================================================

function [v,lam] = rqi(A,x0,ep)
    [m,n] = size(A);
    v = x0;
    lam = v'*A*v;
    I = eye(n);
    kmax = 100;
    for k = 1:kmax
        u = A - lam*I;
        w = u\v;
        v2 = w/norm(w);
        lam2 = v2'*A*v2;
        
        iter = k; 
        
        %convergence
        if norm(A*v2 - lam2*v2) < ep
            break
        end
        v = v2;
        lam = lam2;
        
    end
    fprintf("The code converged at %d iterations to solution\n",iter);

    
end
function [A,absA]=PDEEigenvectors(Q)

global pdetype gamma g;
%
switch(pdetype)
    case(0) % Euler
        % Primitive variables
        rho = Q(1);
        u   = Q(2)/Q(1);
        p   = (gamma-1)*( Q(3)-0.5*rho*u*u );
        % Auxiliary variables
        c   = sqrt(gamma*p/rho);
        H   = (Q(3)+p)/rho;
        v2  = u*u;
        M   = sqrt(v2)/c;
        r2c = rho/2/c;        
        % Jacobian matrix for the Euler equations in conservation form
        A(1,:) = [ 0, 1, 0 ];
        A(2,1) = -u*u+0.5*(gamma-1)*v2;
        A(2,2) = (3.-gamma)*u;
        A(2,3) = gamma-1;
        A(3,1) = -u*(gamma*Q(3)/Q(1)-(gamma-1.)*v2);
        A(3,2) = gamma*Q(3)/Q(1)-0.5*(gamma-1.)*(v2+2*u*u);
        A(3,3) = gamma*u;        
        % Right eigenvector matrix
        R(1,1:3)  = [ 1.,     r2c,         r2c           ];
        R(2,1:3)  = [ u,      r2c*(u+c),   r2c*(u-c)     ];
        R(3,1:3)  = [ 0.5*v2, r2c*(H+c*u), r2c*(H-c*u)   ];        
        % Left eigenvector matrix (inverse of R)
        iR(1,1) = 1.-(gamma-1.)/2.*M*M;
        iR(1,2) =    (gamma-1.)*u/c/c;
        iR(1,3) =   -(gamma-1.)/c/c;
        iR(2,1) = c/rho*(0.5*(gamma-1.)*M*M-u/c);
        iR(2,2) = 1./rho*( 1.-(gamma-1.)*u/c);
        iR(2,3) = (gamma-1.)/rho/c;
        iR(3,1) = c/rho*(0.5*(gamma-1.)*M*M+u/c);
        iR(3,2) = 1./rho*(-1.-(gamma-1.)*u/c);
        iR(3,3) = (gamma-1.)/rho/c;        
         % Eigenvalues
        L(1,1) = u;
        L(2,2) = u+c;
        L(3,3) = u-c;        
        % Compute the absolute value of matrix A
        absA = R*abs(L)*iR;
        
    case(1) % Shallow water equations with scalar transport
        u = Q(2)/Q(1);
        psi = Q(3)/Q(1);
        c = sqrt(g*Q(1)); 
        % Jacobian
        A = [0, 1, 0;
            c^2 - u^2, 2*u, 0;
            -u*psi, psi, u;];
        % Eigenvector matrix
        R = [1, 0, 1;
            u - c, 0, u + c;
            psi, 1, psi;];        
        iR = [(u + c) / c / 2, -1 /(2*c), 0;
            -psi, 0, 1;
            -(u - c) / c / 2,  1/(2*c), 0];
        % Eigenvalues
        L = zeros(3,3);
        L(1,1) = u-c;
        L(2,2) = u;
        L(3,3) = u+c;        
        % Compute the absolute value of matrix A
        absA = R*abs(L)*iR;
        
    case(2) % Shallow water equations with variable bed
        h = Q(1);
        u = Q(2)/h;
        c = sqrt(g*h); 
        % Jacobian
        A = [0,1,0;
            c^2-u^2, 2*u, c^2;
            0, 0, 0];
        % Eigenvector matrix
        R = [1,   1,  1;
             u-c, 0,  u+c;
             0,   (u/c)^2-1, 0 ];     
        iR = 0.5*[(u+c)/c, -1/c, c/(c-u);
                  0,       0,   2*c^2/(u^2-c^2);
                  (c-u)/c,  1/c, c/(c+u)];
        % Eigenvalues
        L = [u-c, 0, 0;
                  0,   0, 0;
                  0,   0, u+c];
        % Compute the absolute value of matrix A
        absA = R*abs(L)*iR;
end
end
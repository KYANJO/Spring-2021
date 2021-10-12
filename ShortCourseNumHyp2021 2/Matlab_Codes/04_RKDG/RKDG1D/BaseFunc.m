function [phi,phixi]=BaseFunc(xi)
global N

switch(N)
    case(0)
        phi(1)=1;               % phi_0
        phixi(1)=0;             % d/dxi phi_0
    case(1)
        phi(1)=1;               % phi_0
        phi(2)=2*xi-1;          % phi_1
        
        phixi(1)=0;             % d/dxi phi_0
        phixi(2)=2;             % d/dxi phi_1
    case(2)
        phi(1)=1;               % phi_0
        phi(2)=2*xi-1;          % phi_1
        phi(3)=1-6*xi+6*xi^2;   % phi_2
        
        phixi(1)=0;             % d/dxi phi_0
        phixi(2)=2;             % d/dxi phi_1
        phixi(3)=-6+12*xi;      % d/dxi phi_2
    case(3)
        phi(1) = 1;
        phi(2) = 2 * xi - 1;
        phi(3) = 6 * xi^2 - 6 * xi + 1;
        phi(4) = 20 * xi^3 - 30 * xi^2 + 12 * xi - 1;
        
        phixi(1) = 0;
        phixi(2) = 2;
        phixi(3) = 12 * xi - 6;
        phixi(4) = 60 * xi^2 - 60 * xi + 12;
    case(4)
        phi(1) = 1;
        phi(2) = 2 * xi - 1;
        phi(3) = 6 * xi^2 - 6 * xi + 1;
        phi(4) = 20 * xi^3 - 30 * xi^2 + 12 * xi - 1;
        phi(5) = 70 * xi^4 - 140 * xi^3 + 90 * xi^2 - 20 * xi + 1;
        
        phixi(1) = 0;
        phixi(2) = 2;
        phixi(3) = 12 * xi - 6;
        phixi(4) = 60 * xi^2 - 60 * xi + 12;
        phixi(5) = 280 * xi^3 - 420 * xi^2 + 180 * xi - 20;                
    case(5)        
        phi(1) = 1;
        phi(2) = 2 * xi - 1;
        phi(3) = 6 * xi^2 - 6 * xi + 1;
        phi(4) = 20 * xi^3 - 30 * xi^2 + 12 * xi - 1;
        phi(5) = 70 * xi^4 - 140 * xi^3 + 90 * xi^2 - 20 * xi + 1;
        phi(6) = 252 * xi^5 - 630 * xi^4 + 560 * xi^3 - 210 * xi^2 + 30 * xi - 1;
        
        phixi(1) = 0;
        phixi(2) = 2;
        phixi(3) = 12 * xi - 6;
        phixi(4) = 60 * xi^2 - 60 * xi + 12;
        phixi(5) = 280 * xi^3 - 420 * xi^2 + 180 * xi - 20;
        phixi(6) = 1260 * xi^4 - 2520 * xi^3 + 1680 * xi^2 - 420 * xi + 30;
        
end
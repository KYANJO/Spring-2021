function [phi, phi_xi, phi_eta] =BaseFunc(xi,eta)
global N
% The orthogonal Dubiner basis on triangles 
switch(N)
    case(0)
        % piecewise constant polynomials of approximation degree N=0 
        % (first order scheme) 
        phi( 1) = 1;
        phi_xi( 1) = 0;
        phi_eta( 1) = 0;
    case(1) 
        % piecewise linear polynomials of approximation degree N=1 
        % (second order scheme) 
        phi( 1) = 1;
        phi( 2) = 2 * xi - 1 + eta;
        phi( 3) = -1 + 3 * eta;
        
        phi_xi( 1) = 0;
        phi_xi( 2) = 2;
        phi_xi( 3) = 0;
        
        phi_eta( 1) = 0;
        phi_eta( 2) = 1;
        phi_eta( 3) = 3;                
    case(2)
        % piecewise quadratic polynomials of approximation degree N=2 
        % (third order scheme)         
        phi( 1) = 1;
        phi( 2) = 2 * xi - 1 + eta;
        phi( 3) = -1 + 3 * eta;
        phi( 4) = eta ^ 2 + 6 * eta * xi + 6 * xi ^ 2 - 2 * eta - 6 * xi + 1;
        phi( 5) = 5 * eta ^ 2 + 10 * eta * xi - 6 * eta - 2 * xi + 1;
        phi( 6) = 10 * eta ^ 2 - 8 * eta + 1;

        phi_xi( 1) = 0;
        phi_xi( 2) = 2;
        phi_xi( 3) = 0;
        phi_xi( 4) = 6 * eta + 12 * xi - 6;
        phi_xi( 5) = 10 * eta - 2;
        phi_xi( 6) = 0;

        phi_eta( 1) = 0;
        phi_eta( 2) = 1;
        phi_eta( 3) = 3;
        phi_eta( 4) = 2 * eta + 6 * xi - 2;
        phi_eta( 5) = 10 * eta + 10 * xi - 6;
        phi_eta( 6) = 20 * eta - 8;
        
    case(3)
        % piecewise cubic polynomials of approximation degree N=3 
        % (fourth order scheme) 
        
        phi( 1) = 1;
        phi( 2) = 2 * xi - 1 + eta;
        phi( 3) = -1 + 3 * eta;
        phi( 4) = eta ^ 2 + 6 * eta * xi + 6 * xi ^ 2 - 2 * eta - 6 * xi + 1;
        phi( 5) = 5 * eta ^ 2 + 10 * eta * xi - 6 * eta - 2 * xi + 1;
        phi( 6) = 10 * eta ^ 2 - 8 * eta + 1;
        phi( 7) = eta ^ 3 + 12 * eta ^ 2 * xi + 30 * eta * xi ^ 2 + 20 * xi ^ 3 - 3 * eta ^ 2 - 24 * eta * xi - 30 * xi ^ 2 + 3 * eta + 12 * xi - 1;
        phi( 8) = 7 * eta ^ 3 + 42 * eta ^ 2 * xi + 42 * eta * xi ^ 2 - 15 * eta ^ 2 - 48 * eta * xi - 6 * xi ^ 2 + 9 * eta + 6 * xi - 1;
        phi( 9) = 21 * eta ^ 3 + 42 * eta ^ 2 * xi - 33 * eta ^ 2 - 24 * eta * xi + 13 * eta + 2 * xi - 1;
        phi(10) = 35 * eta ^ 3 - 45 * eta ^ 2 + 15 * eta - 1;   
        
        phi_xi( 1) = 0;
        phi_xi( 2) = 2;
        phi_xi( 3) = 0;
        phi_xi( 4) = 6 * eta + 12 * xi - 6;
        phi_xi( 5) = 10 * eta - 2;
        phi_xi( 6) = 0;
        phi_xi( 7) = 12 * eta ^ 2 + 60 * eta * xi + 60 * xi ^ 2 - 24 * eta - 60 * xi + 12;
        phi_xi( 8) = 42 * eta ^ 2 + 84 * eta * xi - 48 * eta - 12 * xi + 6;
        phi_xi( 9) = 42 * eta ^ 2 - 24 * eta + 2;
        phi_xi(10) = 0;

        phi_eta( 1) = 0;
        phi_eta( 2) = 1;
        phi_eta( 3) = 3;
        phi_eta( 4) = 2 * eta + 6 * xi - 2;
        phi_eta( 5) = 10 * eta + 10 * xi - 6;
        phi_eta( 6) = 20 * eta - 8;
        phi_eta( 7) = 3 * eta ^ 2 + 24 * eta * xi + 30 * xi ^ 2 - 6 * eta - 24 * xi + 3;
        phi_eta( 8) = 21 * eta ^ 2 + 84 * eta * xi + 42 * xi ^ 2 - 30 * eta - 48 * xi + 9;
        phi_eta( 9) = 63 * eta ^ 2 + 84 * eta * xi - 66 * eta - 24 * xi + 13;
        phi_eta(10) = 105 * eta ^ 2 - 90 * eta + 15;
        
end
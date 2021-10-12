function [theta,theta_xi,theta_tau]=SpaceTimeBaseFunc(xi,tau)

global N;

if(N==0)
    theta(1,1)     = 1;
    theta_xi(1,1)  = 0;
    theta_tau(1,1) = 0;
elseif(N==1)
    theta(1,1) = 1-xi-tau+xi*tau;
    theta(2,1) = xi-xi*tau;
    theta(3,1) = tau-xi*tau;
    theta(4,1) = xi*tau; 
    %                             
    theta_xi(1,1) = -1+tau;
    theta_xi(2,1) = 1-tau;
    theta_xi(3,1) = -tau;
    theta_xi(4,1) = +tau; 
    %
    theta_tau(1,1) = -1+xi;
    theta_tau(2,1) = -xi;
    theta_tau(3,1) = 1-xi;
    theta_tau(4,1) = xi; 
    %
elseif(N==2)
    %
    theta(1,1) = (2*xi-1)*(xi-1)*(2*tau-1)*(tau-1); 
    theta(2,1) = -4*xi*(2*tau-1)*(tau-1)*(xi-1);
    theta(3,1) = xi*(2*tau-1)*(tau-1)*(2*xi-1);
    theta(4,1) = -4*tau*(2*xi-1)*(xi-1)*(tau-1);
    theta(5,1) = 16*tau*xi*(xi-1)*(tau-1);
    theta(6,1) = -4*tau*xi*(2*xi-1)*(tau-1);
    theta(7,1) = tau*(2*xi-1)*(xi-1)*(2*tau-1);
    theta(8,1) = -4*tau*xi*(xi-1)*(2*tau-1);
    theta(9,1) = tau*xi*(2*xi-1)*(2*tau-1);
    %
    theta_xi(1,1) = (2*tau-1)*(tau-1)*(-3+4*xi);
    theta_xi(2,1) = -4*(2*tau-1)*(tau-1)*(2*xi-1);
    theta_xi(3,1) = (2*tau-1)*(tau-1)*(-1+4*xi);
    theta_xi(4,1) = -4*tau*(-3+4*xi)*(tau-1);
    theta_xi(5,1) = 16*tau*(2*xi-1)*(tau-1);
    theta_xi(6,1) = -4*tau*(-1+4*xi)*(tau-1);
    theta_xi(7,1) = tau*(-3+4*xi)*(2*tau-1);
    theta_xi(8,1) = -4*tau*(2*xi-1)*(2*tau-1);
    theta_xi(9,1) = tau*(-1+4*xi)*(2*tau-1);
    %
    theta_tau(1,1) = (2*xi-1)*(xi-1)*(-3+4*tau);
    theta_tau(2,1) = -4*xi*(xi-1)*(-3+4*tau);
    theta_tau(3,1) = xi*(2*xi-1)*(-3+4*tau);
    theta_tau(4,1) = -4*(2*xi-1)*(xi-1)*(2*tau-1);
    theta_tau(5,1) = 16*xi*(xi-1)*(2*tau-1);
    theta_tau(6,1) = -4*xi*(2*xi-1)*(2*tau-1);
    theta_tau(7,1) = (2*xi-1)*(xi-1)*(-1+4*tau);
    theta_tau(8,1) = -4*xi*(xi-1)*(-1+4*tau);
    theta_tau(9,1) = xi*(2*xi-1)*(-1+4*tau);
    %
elseif(N==3)
    %
    theta( 1,1) = (tau-1)*(3*tau-1)*(3*tau-2)*(xi-1)*(3*xi-1)*(3*xi-2)/4;
    theta( 2,1) = -9.D0/4.D0*xi*(tau-1)*(3*tau-1)*(3*tau-2)*(3*xi-2)*(xi-1);
    theta( 3,1) = 9.D0/4.D0*xi*(tau-1)*(3*tau-1)*(3*tau-2)*(3*xi-1)*(xi-1);
    theta( 4,1) = -xi*(tau-1)*(3*tau-1)*(3*tau-2)*(3*xi-1)*(3*xi-2)/4;
    theta( 5,1) = -9.D0/4.D0*tau*(3*tau-2)*(tau-1)*(xi-1)*(3*xi-1)*(3*xi-2);
    theta( 6,1) = 81.D0/4.D0*xi*tau*(3*tau-2)*(tau-1)*(3*xi-2)*(xi-1);
    theta( 7,1) = -81.D0/4.D0*xi*tau*(3*tau-2)*(tau-1)*(3*xi-1)*(xi-1);
    theta( 8,1) = 9.D0/4.D0*xi*tau*(3*tau-2)*(tau-1)*(3*xi-1)*(3*xi-2);
    theta( 9,1) = 9.D0/4.D0*tau*(3*tau-1)*(tau-1)*(xi-1)*(3*xi-1)*(3*xi-2);
    theta(10,1) = -81.D0/4.D0*xi*tau*(3*tau-1)*(tau-1)*(3*xi-2)*(xi-1);
    theta(11,1) = 81.D0/4.D0*xi*tau*(3*tau-1)*(tau-1)*(3*xi-1)*(xi-1);
    theta(12,1) = -9.D0/4.D0*xi*tau*(3*tau-1)*(tau-1)*(3*xi-1)*(3*xi-2);
    theta(13,1) = -tau*(3*tau-1)*(3*tau-2)*(xi-1)*(3*xi-1)*(3*xi-2)/4;
    theta(14,1) = 9.D0/4.D0*xi*tau*(3*tau-1)*(3*tau-2)*(3*xi-2)*(xi-1);
    theta(15,1) = -9.D0/4.D0*xi*tau*(3*tau-1)*(3*tau-2)*(3*xi-1)*(xi-1);
    theta(16,1) = xi*tau*(3*tau-1)*(3*tau-2)*(3*xi-1)*(3*xi-2)/4;
    % 
    theta_xi( 1,1) = (tau-1)*(3*tau-1)*(3*tau-2)*(11-36*xi+27*xi^2)/4;
    theta_xi( 2,1) = -9.D0/4.D0*(tau-1)*(3*tau-1)*(3*tau-2)*(2-10*xi+9*xi^2);
    theta_xi( 3,1) = 9.D0/4.D0*(tau-1)*(3*tau-1)*(3*tau-2)*(1-8*xi+9*xi^2);
    theta_xi( 4,1) = -(tau-1)*(3*tau-1)*(3*tau-2)*(2-18*xi+27*xi^2)/4;
    theta_xi( 5,1) = -9.D0/4.D0*tau*(3*tau-2)*(tau-1)*(11-36*xi+27*xi^2);
    theta_xi( 6,1) = 81.D0/4.D0*tau*(3*tau-2)*(tau-1)*(2-10*xi+9*xi^2);
    theta_xi( 7,1) = -81.D0/4.D0*tau*(3*tau-2)*(tau-1)*(1-8*xi+9*xi^2);
    theta_xi( 8,1) = 9.D0/4.D0*tau*(3*tau-2)*(tau-1)*(2-18*xi+27*xi^2);
    theta_xi( 9,1) = 9.D0/4.D0*tau*(3*tau-1)*(tau-1)*(11-36*xi+27*xi^2);
    theta_xi(10,1) = -81.D0/4.D0*tau*(3*tau-1)*(tau-1)*(2-10*xi+9*xi^2);
    theta_xi(11,1) = 81.D0/4.D0*tau*(3*tau-1)*(tau-1)*(1-8*xi+9*xi^2);
    theta_xi(12,1) = -9.D0/4.D0*tau*(3*tau-1)*(tau-1)*(2-18*xi+27*xi^2);
    theta_xi(13,1) = -tau*(3*tau-1)*(3*tau-2)*(11-36*xi+27*xi^2)/4;
    theta_xi(14,1) = 9.D0/4.D0*tau*(3*tau-1)*(3*tau-2)*(2-10*xi+9*xi^2);
    theta_xi(15,1) = -9.D0/4.D0*tau*(3*tau-1)*(3*tau-2)*(1-8*xi+9*xi^2);
    theta_xi(16,1) = tau*(3*tau-1)*(3*tau-2)*(2-18*xi+27*xi^2)/4;
    % 
    theta_tau( 1,1) = (11-36*tau+27*tau^2)*(xi-1)*(3*xi-1)*(3*xi-2)/4;
    theta_tau( 2,1) = -9.D0/4.D0*xi*(11-36*tau+27*tau^2)*(3*xi-2)*(xi-1);
    theta_tau( 3,1) = 9.D0/4.D0*xi*(11-36*tau+27*tau^2)*(3*xi-1)*(xi-1);
    theta_tau( 4,1) = -xi*(11-36*tau+27*tau^2)*(3*xi-1)*(3*xi-2)/4;
    theta_tau( 5,1) = -9.D0/4.D0*(2-10*tau+9*tau^2)*(xi-1)*(3*xi-1)*(3*xi-2);
    theta_tau( 6,1) = 81.D0/4.D0*xi*(2-10*tau+9*tau^2)*(3*xi-2)*(xi-1);
    theta_tau( 7,1) = -81.D0/4.D0*xi*(2-10*tau+9*tau^2)*(3*xi-1)*(xi-1);
    theta_tau( 8,1) = 9.D0/4.D0*xi*(2-10*tau+9*tau^2)*(3*xi-1)*(3*xi-2);
    theta_tau( 9,1) = 9.D0/4.D0*(1-8*tau+9*tau^2)*(xi-1)*(3*xi-1)*(3*xi-2);
    theta_tau(10,1) = -81.D0/4.D0*xi*(1-8*tau+9*tau^2)*(3*xi-2)*(xi-1);
    theta_tau(11,1) = 81.D0/4.D0*xi*(1-8*tau+9*tau^2)*(3*xi-1)*(xi-1);
    theta_tau(12,1) = -9.D0/4.D0*xi*(1-8*tau+9*tau^2)*(3*xi-1)*(3*xi-2);
    theta_tau(13,1) = -(2-18*tau+27*tau^2)*(xi-1)*(3*xi-1)*(3*xi-2)/4;
    theta_tau(14,1) = 9.D0/4.D0*xi*(2-18*tau+27*tau^2)*(3*xi-2)*(xi-1);
    theta_tau(15,1) = -9.D0/4.D0*xi*(2-18*tau+27*tau^2)*(3*xi-1)*(xi-1);
    theta_tau(16,1) = xi*(2-18*tau+27*tau^2)*(3*xi-1)*(3*xi-2)/4;
end



function [psi,psi_xi]=SpaceBaseFunc(xi)

global N;

if(N==0)
    %
    psi(1,1) = 1;
    %
    psi_xi(1,1) = 0;
    %
elseif(N==1)
    psi(1,1) = 1;
    psi(2,1) = 2*xi-1;
    %
    psi_xi(1,1) = 0;
    psi_xi(2,1) = 2;
elseif(N==2)
    %
    psi(1,1) = 1;
    psi(2,1) = 2*xi-1;
    psi(3,1) = 1-6*xi+6*xi^2;
    %
    psi_xi(2,1) = 2;
    psi_xi(3,1) = -6+12*xi;
    %
elseif(N==3)
    %
    psi( 1,1) = 1;
    psi( 2,1) = 2*xi-1;
    psi( 3,1) = 6*xi^2-6*xi+1;
    psi( 4,1) = 20*xi^3-30*xi^2+12*xi-1;
    %
    psi_xi( 1,1) = 0;
    psi_xi( 2,1) = 2;
    psi_xi( 3,1) = 12*xi-6;
    psi_xi( 4,1) = 60*xi^2-60*xi+12;
    %
elseif(N==4)
    %
    psi( 1,1) = 1;
    psi( 2,1) = 2*xi-1;
    psi( 3,1) = 6*xi^2-6*xi+1;
    psi( 4,1) = 20*xi^3-30*xi^2+12*xi-1;
    psi( 5,1) = 70*xi^4-140*xi^3+90*xi^2-20*xi+1;
    %
    psi_xi( 1,1) = 0;
    psi_xi( 2,1) = 2;
    psi_xi( 3,1) = 12*xi-6;
    psi_xi( 4,1) = 60*xi^2-60*xi+12;
    psi_xi( 5,1) = 4*70*xi^3-3*140*xi^2+2*90*xi-20;
    %
elseif(N==5)
    %
    psi( 1,1) = 1;
    psi( 2,1) = 2*xi-1;
    psi( 3,1) = 6*xi^2-6*xi+1;
    psi( 4,1) = 20*xi^3-30*xi^2+12*xi-1;
    psi( 5,1) = 70*xi^4-140*xi^3+90*xi^2-20*xi+1;
    psi( 6,1) = 252*xi^5-630*xi^4+560*xi^3-210*xi^2+30*xi-1;
    %
    psi_xi( 1,1) = 0;
    psi_xi( 2,1) = 2;
    psi_xi( 3,1) = 12*xi-6;
    psi_xi( 4,1) = 60*xi^2-60*xi+12;
    psi_xi( 5,1) = 4*70*xi^3-3*140*xi^2+2*90*xi-20;
    psi_xi( 6,1) = 5*252*xi^4-4*630*xi^3+3*560*xi^2-2*210*xi+30;
    %
end


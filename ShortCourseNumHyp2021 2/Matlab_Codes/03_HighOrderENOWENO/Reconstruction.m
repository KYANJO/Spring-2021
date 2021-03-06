function w = Reconstruction(u)
% WENO reconstruction, polynomial w is obtained from the cell averages
global N CWENO
sz = size(u);
IMAX = sz(2);
epsilon = 1e-14;
r = 8;
lambdaC = 1e5;
lambdaL = 1;
lambdaR = 1;

if (CWENO == 1)
    lsum=lambdaC+lambdaL+lambdaR;
    lambdaC = lambdaC/lsum; 
    lambdaL = lambdaL/lsum; 
    lambdaR = lambdaR/lsum; 
end

if (N==0)
    w = u;
end

if(N==2 && CWENO==0)
    % polynomial ENO-type WENO reconstruction
    for i=1:IMAX
        im1 = idx(i-1);
        im2 = idx(i-2); 
        ip1 = idx(i+1);
        ip2 = idx(i+2);
        % Compute reconstruction polynomials for each stencil
        % Central stencil
        a0C = u(:,i); 
        a1C = -1/4*u(:,im1)+1/4*u(:,ip1);
        a2C = -1/6*u(:,i)+1/12*u(:,im1)+1/12*u(:,ip1);
        % Left stencil
        a0L = u(:,i);
        a1L =  3/4*u(:,i)+ 1/4*u(:,im2)-u(:,im1);
        a2L = 1/12*u(:,i)+1/12*u(:,im2)-1/6*u(:,im1);
        % Right stencil
        a0R = u(:,i);
        a1R = -3/4*u(:,i)-1/4*u(:,ip2)+u(:,ip1);
        a2R = 1/12*u(:,i)+1/12*u(:,ip2)-1/6*u(:,ip1);
        % Compute the oscillation indicators
        OIC = 156*a2C.^2+4*a1C.^2;
        OIL = 156*a2L.^2+4*a1L.^2;
        OIR = 156*a2R.^2+4*a1R.^2;
        % Compute the non-normalized nonlinear weights
        omegaCt = lambdaC./(epsilon+OIC).^r;
        omegaLt = lambdaL./(epsilon+OIL).^r;
        omegaRt = lambdaR./(epsilon+OIR).^r;
        % Compute the sum of the weights
        omegasum = omegaCt + omegaLt + omegaRt;
        % Normalize with the sum
        omegaC  = omegaCt./omegasum;
        omegaL  = omegaLt./omegasum;
        omegaR  = omegaRt./omegasum;
         % Compute the final coefficients a,b,c:
        a0 = omegaC.*a0C + omegaL.*a0L + omegaR.*a0R;
        a1 = omegaC.*a1C + omegaL.*a1L + omegaR.*a1R;
        a2 = omegaC.*a2C + omegaL.*a2L + omegaR.*a2R;
        % Put the coefficients into a matrix 
        w(1,:,i) = a0;
        w(2,:,i) = a1;
        w(3,:,i) = a2;
    end  
end

if(N==2 && CWENO==1)
    % CWENO reconstruction (optimal stencil, arbitrary linear weights) 
    for i=1:IMAX
        im1 = idx(i-1);
        ip1 = idx(i+1);
        % Compute reconstruction polynomials for each stencil
        % Popt polynomial (degree 2) 
        a0C = u(:,i); 
        a1C = -1/4*u(:,im1)+1/4*u(:,ip1);
        a2C = -1/6*u(:,i)+1/12*u(:,im1)+1/12*u(:,ip1);
        % Left stencil (degree 1)
        a0L = u(:,i);
        a1L = 1/2*u(:,i)-1/2*u(:,im1);
        a2L = 0; 
        % Right stencil
        % Right stencil (degree 1) 
        a0R = u(:,i);
        a1R = 1/2*u(:,ip1)-1/2*u(:,i);
        a2R = 0; 
        % P0 polynomial 
        lsum = lambdaC+lambdaL+lambdaR; 
        a00 = 1/lambdaC*(lsum*a0C - lambdaL*a0L - lambdaR*a0R); 
        a10 = 1/lambdaC*(lsum*a1C - lambdaL*a1L - lambdaR*a1R); 
        a20 = 1/lambdaC*(lsum*a2C - lambdaL*a2L - lambdaR*a2R); 
         % Compute the oscillation indicators
        OI0 = 156*a20.^2+4*a10.^2;
        OIL = 156*a2L.^2+4*a1L.^2;
        OIR = 156*a2R.^2+4*a1R.^2;
        % Compute the non-normalized nonlinear weights
        omega0t = lambdaC./(epsilon+OI0).^r;
        omegaLt = lambdaL./(epsilon+OIL).^r;
        omegaRt = lambdaR./(epsilon+OIR).^r;
        % Compute the sum of the weights
        omegasum = omega0t + omegaLt + omegaRt;
        % Normalize with the sum
        omega0  = omega0t./omegasum;
        omegaL  = omegaLt./omegasum;
        omegaR  = omegaRt./omegasum;
         % Compute the final coefficients a,b,c:
        a0 = omega0.*a0C + omegaL.*a0L + omegaR.*a0R;
        a1 = omega0.*a1C + omegaL.*a1L + omegaR.*a1R;
        a2 = omega0.*a2C + omegaL.*a2L + omegaR.*a2R;
        % Put the coefficients into a matrix 
        w(1,:,i) = a0;
        w(2,:,i) = a1;
        w(3,:,i) = a2;
    end  
end

if(N==1)
    for i=[1:IMAX]
        im1 = idx(i-1);
        ip1 = idx(i+1);
        % Left stencil
        a0L = u(:,i);
        a1L = 1/2*u(:,i)-1/2*u(:,im1);
        % Right stencil
        a0R = u(:,i);
        a1R = 1/2*u(:,ip1)-1/2*u(:,i);
        % Minmod
        a0 = a0L;
        a1 = minmod(a1L,a1R); 
        w(1,:,i) = a0;   
        w(2,:,i) = a1;
    end
end

if(N==3 && CWENO==0)
    % polynomial ENO-type WENO reconstruction 
      for i=[1:IMAX]
        im1 = idx(i-1);
        im2 = idx(i-2); 
        im3 = idx(i-3);  
        ip1 = idx(i+1);
        ip2 = idx(i+2);
        ip3 = idx(i+3);
        % Compute the reconstruction polynomials for each stencil
        % Central stencil 1
        a0C1 = u(:,i); 
        a1C1 = -21/40*u(:,im1)+11/40*u(:,i)+11/120*u(:,im2)+19/120*u(:,ip1);
        a2C1 = -1/6*u(:,i)+1/12*u(:,im1)+1/12*u(:,ip1);
        a3C1 = 1/40*u(:,im1)-1/40*u(:,i)-1/120*u(:,im2)+1/120*u(:,ip1); 
        % Central stencil 2
        a0C2 = u(:,i); 
        a1C2 = -19/120*u(:,im1)-11/40*u(:,i)-11/120*u(:,ip2)+21/40*u(:,ip1);
        a2C2 = -1/6*u(:,i)+1/12*u(:,im1)+1/12*u(:,ip1);
        a3C2 = -1/120*u(:,im1)+1/40*u(:,i)+1/120*u(:,ip2)-1/40*u(:,ip1); 
        % Left stencil
        a0L = u(:,i);
        a1L = -19/120*u(:,im3)+109/120*u(:,i)+29/40*u(:,im2)-59/40*u(:,im1);
        a2L = -1/12*u(:,im3)+1/6*u(:,i)+1/3*u(:,im2)-5/12*u(:,im1);
        a3L = -1/120*u(:,im3)+1/120*u(:,i)+1/40*u(:,im2)-1/40*u(:,im1); 
        % Right stencil
        a0R = u(:,i);
        a1R = 19/120*u(:,ip3)-109/120*u(:,i)-29/40*u(:,ip2)+59/40*u(:,ip1);
        a2R = -1/12*u(:,ip3)+1/6*u(:,i)+1/3*u(:,ip2)-5/12*u(:,ip1);
        a3R = 1/120*u(:,ip3)-1/120*u(:,i)-1/40*u(:,ip2)+1/40*u(:,ip1); 
        % Compute the oscillation indicators
        OIC1 = 1224*a3C1.^2+8*a3C1.*a1C1+156*a2C1.^2+4*a1C1.^2;
        OIC2 = 1224*a3C2.^2+8*a3C2.*a1C2+156*a2C2.^2+4*a1C2.^2;
        OIL  = 1224*a3L.^2 +8*a3L.*a1L  +156*a2L.^2 +4*a1L.^2;
        OIR  = 1224*a3R.^2 +8*a3R.*a1R  +156*a2R.^2 +4*a1R.^2; 
        % Compute the non-normalized nonlinear weights
        omegaC1t = lambdaC./(epsilon+OIC1).^r;
        omegaC2t = lambdaC./(epsilon+OIC2).^r;
        omegaLt = lambdaL./(epsilon+OIL).^r;
        omegaRt = lambdaR./(epsilon+OIR).^r;
        % Compute the sum of the weights
        omegasum = omegaC1t + omegaC2t + omegaLt + omegaRt;
        % Normalize with the sum
        omegaC1  = omegaC1t./omegasum;
        omegaC2  = omegaC2t./omegasum;
        omegaL   = omegaLt./omegasum;
        omegaR   = omegaRt./omegasum;
        % Compute the final coefficients a,b,c:
        a0 = omegaC1.*a0C1 + omegaC2.*a0C2 + omegaL.*a0L + omegaR.*a0R;
        a1 = omegaC1.*a1C1 + omegaC2.*a1C2 + omegaL.*a1L + omegaR.*a1R;
        a2 = omegaC1.*a2C1 + omegaC2.*a2C2 + omegaL.*a2L + omegaR.*a2R;
        a3 = omegaC1.*a3C1 + omegaC2.*a3C2 + omegaL.*a3L + omegaR.*a3R;
        % Put the coefficients into a matrix 
        w(1,:,i) = a0;
        w(2,:,i) = a1;
        w(3,:,i) = a2;
        w(4,:,i) = a3;
      end
end

if(N==3 && CWENO==1)
    % CWENO reconstruction (optimal stencil, arbitrary linear weights) 
      for i=[1:IMAX]
        im1 = idx(i-1);
        im2 = idx(i-2); 
        ip1 = idx(i+1);
        ip2 = idx(i+2);
        % Compute the reconstruction polynomials for each stencil
        % Central stencil (Popt, degree 3) 
        a0C1 = u(:,i); 
        a1C1 = -21/40*u(:,im1)+11/40*u(:,i)+11/120*u(:,im2)+19/120*u(:,ip1);
        a2C1 = -1/6*u(:,i)+1/12*u(:,im1)+1/12*u(:,ip1);
        a3C1 = 1/40*u(:,im1)-1/40*u(:,i)-1/120*u(:,im2)+1/120*u(:,ip1); 
        % Central stencil 2
        a0C2 = u(:,i); 
        a1C2 = -19/120*u(:,im1)-11/40*u(:,i)-11/120*u(:,ip2)+21/40*u(:,ip1);
        a2C2 = -1/6*u(:,i)+1/12*u(:,im1)+1/12*u(:,ip1);
        a3C2 = -1/120*u(:,im1)+1/40*u(:,i)+1/120*u(:,ip2)-1/40*u(:,ip1); 
        % take average of the two 
        a0C = 0.5*(a0C1+a0C2); 
        a1C = 0.5*(a1C1+a1C2); 
        a2C = 0.5*(a2C1+a2C2); 
        a3C = 0.5*(a3C1+a3C2); 
        % Left stencil (degree 1)
        a0L = u(:,i);
        a1L = 1/2*u(:,i)-1/2*u(:,im1);
        a2L = 0; 
        a3L = 0; 
        % Right stencil (degree 1) 
        a0R = u(:,i);
        a1R = 1/2*u(:,ip1)-1/2*u(:,i);
        a2R = 0; 
        a3R = 0; 
        % P0 polynomial 
        lsum = lambdaC+lambdaL+lambdaR; 
        a00 = 1/lambdaC*(lsum*a0C - lambdaL*a0L - lambdaR*a0R); 
        a10 = 1/lambdaC*(lsum*a1C - lambdaL*a1L - lambdaR*a1R); 
        a20 = 1/lambdaC*(lsum*a2C - lambdaL*a2L - lambdaR*a2R);         
        a30 = 1/lambdaC*(lsum*a3C - lambdaL*a3L - lambdaR*a3R);         
        % Compute the oscillation indicators
        OI0  = 1224*a30.^2 +8*a30.*a10  +156*a20.^2 +4*a10.^2;
        OIL  = 1224*a3L.^2 +8*a3L.*a1L  +156*a2L.^2 +4*a1L.^2;
        OIR  = 1224*a3R.^2 +8*a3R.*a1R  +156*a2R.^2 +4*a1R.^2; 
        % Compute the non-normalized nonlinear weights
        omega0t = lambdaC./(epsilon+OI0).^r;
        omegaLt = lambdaL./(epsilon+OIL).^r;
        omegaRt = lambdaR./(epsilon+OIR).^r;
        % Compute the sum of the weights
        omegasum = omega0t + omegaLt + omegaRt;
        % Normalize with the sum
        omega0  = omega0t./omegasum;
        omegaL  = omegaLt./omegasum;
        omegaR  = omegaRt./omegasum;
        % Compute the final coefficients 
        a0 = omega0.*a0C + omegaL.*a0L + omegaR.*a0R;
        a1 = omega0.*a1C + omegaL.*a1L + omegaR.*a1R;
        a2 = omega0.*a2C + omegaL.*a2L + omegaR.*a2R;
        a3 = omega0.*a3C + omegaL.*a3L + omegaR.*a3R;
        % Put the coefficients into a matrix 
        w(1,:,i) = a0;
        w(2,:,i) = a1;
        w(3,:,i) = a2;
        w(4,:,i) = a3;
      end
end

end
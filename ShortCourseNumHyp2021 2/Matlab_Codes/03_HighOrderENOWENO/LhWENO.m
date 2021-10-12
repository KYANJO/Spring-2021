function dqdt = LhWENO(q)
% Optimal WENO reconstruction with optimal weights
global N IMAX dx nVar flux
epsilon = 1.e-7;

switch(N)
    case(2)
        % second order TVD scheme with minmod limiter 
        for i=1:IMAX
            im1 = idx(i-1);
            ip1 = idx(i+1);
            slopeR = (q(:,i)-q(:,im1))/dx; 
            slopeL = (q(:,ip1)-q(:,i))/dx; 
            slope  = minmod(slopeL,slopeR); 
            wL(:,i) = q(:,i) - 0.5*dx*slope; 
            wR(:,i) = q(:,i) + 0.5*dx*slope; 
        end
    case(3)
        % third order classical WENO scheme
        r = 2;
        % optimal weights for the point-value left
        lambda0L = 2/3;
        lambda0R = 1/3;
        % optimal weights for the point-value right
        lambda1L = 1/3;
        lambda1R = 2/3; 
        for i=1:IMAX
            im1 = idx(i-1);
            ip1 = idx(i+1);
            % construct the two lower order polynomials of degree 1
            % left stencil
            l0 = 1/2*q(:,im1)+1/2*q(:,i);
            l1 = q(:,i)-q(:,im1);
            % right stencil
            r0 = 3/2*q(:,i)-1/2*q(:,ip1);
            r1 = -q(:,i)+q(:,ip1);
            % Compute the oscillation indicators
            OIL = l1(:).^2;
            OIR = r1(:).^2;
            % ---------------------------------------------
            % reconstruction on the left boundary
            % ---------------------------------------------
            % Compute the non-normalized nonlinear weights
            omegaLt = lambda0L./(epsilon+OIL).^r;
            omegaRt = lambda0R./(epsilon+OIR).^r;
            % Compute the sum of the weights
            omegasum = omegaLt + omegaRt;
            % Normalize with the sum
            omegaL  = omegaLt./omegasum;
            omegaR  = omegaRt./omegasum;
            wL(:,i) = omegaL.*l0 + omegaR.*r0;
            % ---------------------------------------------
            % reconstruction on the right boundary
            % ---------------------------------------------
            % Compute the non-normalized nonlinear weights
            omegaLt = lambda1L./(epsilon+OIL).^r;
            omegaRt = lambda1R./(epsilon+OIR).^r;
            % Compute the sum of the weights
            omegasum = omegaLt + omegaRt;
            % Normalize with the sum
            omegaL  = omegaLt./omegasum;
            omegaR  = omegaRt./omegasum;
            wR(:,i) = omegaL.*(l0+l1) + omegaR.*(r0+r1);
        end        
    case(5)
        % fifth order classical WENO scheme
        r = 4;
        % optimal weights for the point-value left
        lambda0L = 3/10;
        lambda0C = 6/10;
        lambda0R = 1/10;
        % optimal weights for the point-value right
        lambda1L = 1/10;
        lambda1C = 6/10;
        lambda1R = 3/10;
        for i = 1:IMAX
            im2 = idx(i-2);
            im1 = idx(i-1);
            ip1 = idx(i+1);
            ip2 = idx(i+2);
            % Construct the three lower order polynomials f degree 2
            % left stencil
            l0 = 1/3*q(:,i)+5/6*q(:,im1)-1/6*q(:,im2);
            l1 = q(:,i) -q(:,im1);
            l2 = 1/2*q(:,i)-q(:,im1)+1/2*q(:,im2);
            % cetral stencil
            c0 = 5/6*q(:,i) +1/3*q(:,im1) -1/6*q(:,ip1);
            c1 = q(:,i) -q(:,im1);
            c2 = -q(:,i) +1/2*q(:,im1)+1/2*q(:,ip1);
            % right stencil
            r0 = 1/3*q(:,ip2) +11/6*q(:,i)-7/6*q(:,ip1);
            r1 = -q(:,ip2)-2*q(:,i)+3*q(:,ip1);
            r2 = 1/2*q(:,ip2) +1/2*q(:,i)-q(:,ip1);
            % Oscillator indicators
            OIL = (16/3)*l2(:).^2 + 2*l1(:).*l2(:) + l1(:).^2;
            OIC = (16/3)*c2(:).^2 + 2*c1(:).*c2(:) + c1(:).^2;
            OIR = (16/3)*r2(:).^2 + 2*r1(:).*r2(:) + r1(:).^2;
            % Reconstruction on the left boundary
            % Compute non-normalized nonlinear weights
            omegaCt = lambda0C./(epsilon+OIC).^r;
            omegaLt = lambda0L./(epsilon+OIL).^r;
            omegaRt = lambda0R./(epsilon+OIR).^r;
            % Compute the sum of weights
            omegasum = omegaCt + omegaLt + omegaRt;
            % Normalice with sum
            omegaC = omegaCt./omegasum;
            omegaL = omegaLt./omegasum;
            omegaR = omegaRt./omegasum;
            wL(:,i) = omegaL.*l0+omegaC.*c0 +omegaR.*r0;
            % Right boundary
            % Compute non-normalized nonlinear weights
            omegaCt = lambda1C./(epsilon+OIC).^r;
            omegaLt = lambda1L./(epsilon+OIL).^r;
            omegaRt = lambda1R./(epsilon+OIR).^r;
            % Compute the sum of weights
            omegasum = omegaCt + omegaLt + omegaRt;
            % Normalice with sum
            omegaC = omegaCt./omegasum;
            omegaL = omegaLt./omegasum;
            omegaR = omegaRt./omegasum;
            wR(:,i) = omegaL.*(l0+l1+l2)+omegaC.*(c0+c1+c2) +omegaR.*(r0+r1+r2);
        end
end

dqdt = zeros(nVar,IMAX);
for i = 2:IMAX-1
    ip1 = idx(i+1);
    im1 = idx(i-1);
    if (flux == 0) % Rusanov
        % Define maximum signal speed smax
        smax = max( max(abs(PDEEigenvalues(wR(:,i)))), max(abs(PDEEigenvalues(wL(:,ip1)))));
        fp = 0.5*( PDEFlux(wL(:,ip1))' + PDEFlux(wR(:,i))' ) ...
            - 0.5*smax*( wL(:,ip1)-wR(:,i) );
        smax = max( max(abs(PDEEigenvalues(wR(:,im1)))), max(abs(PDEEigenvalues(wL(:,i)))));
        fm = 0.5*( PDEFlux(wL(:,i))' + PDEFlux(wR(:,im1))' ) ...
            - 0.5*smax*( wL(:,i)-wR(:,im1) );
    end
    dqdt(:,i) = -1/dx*(fp-fm);
end
end
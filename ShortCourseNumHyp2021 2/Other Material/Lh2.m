function dqdt=Lh2(q)
global IMAX dx a 
dqdt = zeros(IMAX,1);
% slope reconstruction 
for i=1:IMAX
    im1 = idx(i-1);
    ip1 = idx(i+1);
    % calculate the limited slope using the minmod limiter 
    slope = minmod( (q(ip1)-q(i))/dx, (q(i)-q(im1))/dx );
    % calculate the boundary-extrapolated values 
    wR(i) = q(i) + 0.5*dx*slope; 
    wL(i) = q(i) - 0.5*dx*slope; 
end
% flux calculation using the upwind scheme 
for i=1:IMAX
    im1 = idx(i-1);
    ip1 = idx(i+1);
    fp = 0.5*a*(wR(i)+wL(ip1))-0.5*abs(a)*(wL(ip1)-wR(i)); 
    fm = 0.5*a*(wR(im1)+wL(i))-0.5*abs(a)*(wL(i)-wR(im1)); 
    dqdt(i) = - 1/dx*( fp - fm );
end

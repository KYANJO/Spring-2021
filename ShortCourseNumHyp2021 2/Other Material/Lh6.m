function dqdt=Lh6(q)
global IMAX dx a 
dqdt = zeros(IMAX,1);
for i=1:IMAX
    % right boundary extrapolated value of cell i 
    im1 = idx(i-1);
    im2 = idx(i-2);
    im3 = idx(i-3);
    i0  = idx(i);
    ip1 = idx(i+1);
    ip2 = idx(i+2);
    qp  = -(1/30)*q(ip2)-(1/60)*q(im3)+(7/60)*q(im2)-(23/60)*q(im1)+(19/20)*q(i0)+(11/30)*q(ip1);   
    fp = a*qp;
    % right boundary extrapolated value of cell i-1 
    im1 = idx(i-1-1);
    im2 = idx(i-2-1);
    im3 = idx(i-3-1);
    i0  = idx(i-1);
    ip1 = idx(i+1-1);
    ip2 = idx(i+2-1);
    qm  = -(1/30)*q(ip2)-(1/60)*q(im3)+(7/60)*q(im2)-(23/60)*q(im1)+(19/20)*q(i0)+(11/30)*q(ip1);   
    fm = a*qm;
    dqdt(i) = - 1/dx*( fp - fm );
end

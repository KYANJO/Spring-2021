function [flux,D]=RusanovFluxWB(QL,QR)
% Compute Rusanov flux well-balanced
LL = lambda(QL);
LR = lambda(QR);

smax = max(max(abs(LL)),max(abs(LR)));

IWB = [1, 0, 1;
    0, 1, 0;
    0, 0, 0];       % Well-balanced identity matrix

flux = 0.5*(f(QL)+f(QR)) -0.5*smax*IWB*(QR-QL);

Qav = 0.5*(QR+QL);  % Average
Btilde = B(Qav);    % Roe matrix for B (path integral)

D = 0.5*Btilde*(QR-QL); % Path conservative jump term
end
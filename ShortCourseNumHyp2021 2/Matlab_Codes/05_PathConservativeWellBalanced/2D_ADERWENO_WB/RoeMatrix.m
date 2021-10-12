function BRoe = RoeMatrix(QL,QR)
% Compute Roe matrix
global nVar
sGP = [0.5-sqrt(15)/10,0.5,0.5+sqrt(15)/10];
wGP = [5/18,8/18,5/18];
BRoe = zeros(nVar,nVar);
for i=1:length(sGP)
    QGP = QL + sGP(i)*(QR-QL); % Linear segment path in Gauss point i
    BGP = PDEMatrixB(QGP); % Matrix B containing the non-conservative terms
    BRoe = BRoe + wGP(i)*BGP;
end
end
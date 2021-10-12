function flux=Rusanovp(QL,QR)
% Rusanov flux, not well balanced
LL=lambda(QL); 
LR=lambda(QR);
smax = max( max(abs(LL)), max(abs(LR)) ); 
flux = 0.5*( fpr(QL) + fpr(QR) ) - 0.5*smax*( QR - QL ); 
end
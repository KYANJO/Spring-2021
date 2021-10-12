% condizione iniziale 
function q=h(x) 
qL = 2; 
qR = 1; 
% problema di Riemann 
% if(x<0) 
%     q = qL; 
% else
%     q = qR; 
% end
% Gaussiana 
q = exp(-0.5*(x+0.0)^2/0.05^2); 
%lambda = 1/5;    % lunghezza dell'onda 
%q = sin(2*pi/lambda*x); 

%q = 1 + 0.1*x^4/24; 
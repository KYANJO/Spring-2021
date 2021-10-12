function q = h(x)
% Condizione iniziale
% Test1: RP
% qL = 2; 
% qR = 1;
% if (x<0)
%     q = qL;
% else
%     q = qR;
% end

% Test2: 
 q = 1+ exp(-0.5*x^2/0.02^2);

end
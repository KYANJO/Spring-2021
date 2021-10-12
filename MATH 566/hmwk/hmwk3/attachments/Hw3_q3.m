clear all
close all
clc

J = 0:520;

x = 2.^(-52+J./10);

ya = log(x+1)./x;

figure
semilogx(x,ya,'r','linewidth',2)
%set(gca,'FontSize',13)
hold on
z = 1 + x;
y = log(z)./(z-1);

semilogx(x,y,'b--','linewidth',2)
xlabel('x')
ylabel('f(x)')
title('Plot of function f(x) = log(x+1)/x')
legend('Direct computation', 'Second method')
set(gca,'FontSize',13)
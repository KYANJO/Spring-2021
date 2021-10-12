clear all
close all

%3a). Compute the condition number for values very close to zero.
%condition number
C = @(x) abs((x - (x+1)*log(x + 1))/((x+1)*log(x + 1)));

n = 10;

xl = linspace(-0.05,-0.0000001,n);
xr = linspace(0.0000001,0.05,n);

cl = zeros(n,1);
cr = zeros(n,1);
for i= 1:n
   cl(i) = C(xl(i));
    cr(i) = C(xr(i));
end

%What does the condition number tell you about the stability of evaluating
%f(x) near zero.
Table = table(xl',cl,xr',cr, 'VariableNames',{'x<0','C(x<0)','x>0','C(x>0)'})
fprintf('The value of x near zero for the condition number, are all small meaning the function is well condition and stable,\n since a small input to the condition number yeild a small output, as observed from the table for different values of x\n near zero\n');

fprintf('3b.\n');
%Evaluate th function f(x) = log(x +1)/x using the expression as given for x
f = @(x) (log(x+1))./x;

j = [0:520]';
xj = 2.^(-52 + j./10);
fj = f(xj);

%plot of f
semilogx(xj,fj);
hold on
title('f & z against x');
xlabel('x'); ylabel('f');
fprintf('The algorithm looks to be unstable near x = 0, according to the distortion of the curve observed \n near that point \n');

fprintf('3c.\n');
%Now evaluate f(x) at the same xj values as part (b)
z = 1 + xj;
y = log(z)./(z-1);
semilogx(xj,y);
legend('f(x)','z');

fprintf('Near x = 0, their is no noise the curve is stable, but in part (b), there is alot of noise in the region \n'); 
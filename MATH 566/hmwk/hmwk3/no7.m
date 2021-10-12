clear all;
close all;

%Skeel condition number (CN).

% identity matrix
I = eye(10);
%Permutation Matrix P
P = [I(:,4) I(:,7) I(:,8) I(:,5) I(:,2) I(:,9) I(:,10)...
    I(:,3) I(:,6) I(:,1)];

%Verify that this is true for both the standard & Skeel 
%standard CN
Cs = cond(P)
fprintf('Hence the standard condition number for P is 1\n');

%Skeel CN
Sc = norm((abs(inv(P))*abs(P)),2)
fprintf('Hence the skeel condition number for P is 1\n');

%Scaling the third column of P
P(:,3) = (10^(-10))*P(:,3);
%standard CN
Cs = cond(P)
fprintf('Hence the standard condition number for P after scaling 1s: %e \n',Cs);

%Skeel CN
Sc = norm((abs(inv(P))*abs(P)),2)
fprintf('Hence the skeel condition number for P after scaling is: 1\n');

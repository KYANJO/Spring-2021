clear all
close all
clc

I = eye(10);
P = [I(:,4),I(:,7),I(:,8),I(:,5),I(:,2),I(:,9),I(:,10),I(:,3),I(:,6),I(:,1)];

% Verification of the condition number of permutation matrix P
fprintf('Condition number using standard condition:\n\n');
cdstd = cond(P);
disp(cdstd)
fprintf('Condition number using Skeel condition:\n\n');
KCR = norm(abs(inv(P))*abs(P));
disp(KCR)

% Scaling the third column of P by 10^(−10)
P(:,3) = 10^(-10)*P(:,3);

fprintf('Condition number using standard condition after scaling:\n\n');
cdstdn = cond(P);
disp(cdstdn)
fprintf('Condition number using Skeel condition after scaling:\n\n');
KCRn = norm(abs(inv(P))*abs(P));
disp(KCRn)

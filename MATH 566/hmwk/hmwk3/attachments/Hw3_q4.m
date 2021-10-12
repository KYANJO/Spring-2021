clear all
close all
clc

A = 1./hankel(2:6,6:10);
b = [0.882 0.744 0.618 0.521 0.447]';
x = A\b;

cd = cond(A);
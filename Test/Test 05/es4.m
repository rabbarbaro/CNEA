clear
clc
close all

A = [8 -2 -2
    -2 6 -1
    -2 -1 9];
b = ones(3, 1);
x0 = b;

x = pcg(A, b, 1e-12, 2, eye(3), eye(3), x0)
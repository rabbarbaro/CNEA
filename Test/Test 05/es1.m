clear
clc
close all

A = [4 2 1
    2 4 1
    1 1 7];
b = 2*ones(3, 1);

x0 = b;

% gradiente / -r
r_0 = b - A*x0;
alpha_0 = (r_0' * r_0) / (r_0' * A * r_0)

x1 = x0 + alpha_0 * r_0
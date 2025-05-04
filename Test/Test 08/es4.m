clear
clc
close all

x_t = [997, 1975, 2153, 2722, 2431, 2546, 1782, 1040, 1670];

j = 5;
n = 3;

x_tj = x_t(j-(-n:n));

x_tilde_t = 1/(2*n + 1) * sum(x_tj)
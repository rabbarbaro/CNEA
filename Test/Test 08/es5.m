clear
clc
close all

f = @(x) exp(x) .* sin(pi*x);

a = -1;
b = 1;

x_q = linspace(a, b, 1000);

x_nod = linspace(a, b, 11);
f_nod = f(x_nod);

pz = interp1(x_nod, f_nod, 0.7)
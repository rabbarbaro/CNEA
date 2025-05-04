clear
clc
close all

x = linspace(0, 1, 150);
rng(1);
y = 3 * x.^2 + 0.3 * sin(100 * pi * x) + 0.3 * randn(1, 150);

x_q = linspace(0, 1, 1000);

p = polyfit(x, y, 2);
p_disp = polyval(p, x_q);

plot(x, y, x_q, p_disp)

pval = polyval(p, 1.5)
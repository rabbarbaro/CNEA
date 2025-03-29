clear
clc
close all

% definisco la funzione come due distinte
f1 = @(x) sin(x);
f2 = @(x) (x.^2 + x)./6;

% definisco due intervalli distinti
x1 = -pi:0.01:0;
x2 = 0:0.01:pi;

plot(x1, f1(x1), x2, f2(x2))
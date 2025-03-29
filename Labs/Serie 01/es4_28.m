clear
clc
close all

% definisco la funzione
f = @(x) (x.^2 + 3)./(x - 1);

% plotto
x = -5:0.01:5;
plot(x, f(x))
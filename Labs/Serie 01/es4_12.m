clear
clc
close all

% definisco la funzione a tratti moltiplicando i vari tratti per il vettore
% logico in uscita dalla verifica dei tratti
f = @(x) -sqrt(x.^2 - x).*(x < 0) + (-x.^2 + 2*x).*exp(-x).*(x >= 0);
x = -5:.01:5;
plot(x, f(x));
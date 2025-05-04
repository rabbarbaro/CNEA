clear
clc
close all

a = -1;
b = 1;

f = @(x) exp(-x) - beta * x + gamma;
df = @(x) -exp(-x) - beta;
d2f = @(x) exp(-x);

% dalla teoria l'errore per il punto medio composito
% err_PMcomp = (b-a)/24 * H^2 * d2f(xi);

tol = 1e-2;
d2f_max = exp(1);
H_max = sqrt(tol * 24 / ((b-a) * d2f_max));
N_min = (b-a) / H_max
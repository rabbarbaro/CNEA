clear
clc
close all

f = @(x) cos(pi * x);

a = -1;
b = 1;

% dalla teoria la stima a priori dell'errore di interpolazione con lagrange su n+1 nodi equispaziati è:
%   e_n = 1/(4 * (n+1)) * ((b-a)/n)^(n+1) * max|d(n+1)f(x)|

% la derivata di f n-esima è pi^n * funzione trigonometrica limitata a -1/1

% quindi:
%   e_n = 1/(4 * (n+1)) * (2/n)^(n+1) * pi^n+1
%   e_n = 2^(n-1) / (n+1) * (pi/n)^(n+1)
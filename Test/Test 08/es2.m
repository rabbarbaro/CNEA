clear
clc
close all

f = @(x) cos(pi * x);

a = -1;
b = 1;

% dalla teoria la stima a priori dell'errore di interpolazione con lagrange su n+1 nodi equispaziati è:
%   - H = (b-a) / n
%   - E_H <= H^2 / 8 * max|f''(x)|, x in [a,b]

d2f = @(x) -pi^2 * cos(pi * x); % derivata seconda di f
max_d2f = max(abs(d2f(linspace(a, b, 100)))); % massimo della derivata seconda di f in [a,b] (è pi^2)
H = (b - a) / 2; % passo di interpolazione
E_H = H^2 / 8 * max_d2f; % errore di interpolazione
disp(['Errore di interpolazione: ', num2str(E_H)]);
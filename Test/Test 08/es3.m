clear
clc
close all

n = 3;
a = 0;
b = 1;

% calcolo la base di polinomi caratteristici di lagrange di grado n = 3

% creo un vettore con i parametri t_k per k = 0, ..., n
% di fatto sono gli n+1 nodi di Chebyshev nell'intervallo di
% riferimento [-1, 1]
tv = -cos(pi * (0:n)/n);
% mappo i nodi dall'intervallo di riferimento ad [a, b]
x_nod_che = (b - a)/2 * tv + (a + b)/2;
% definisco i punti x_v per i nodi di Chebyshev
x_v = x_nod_che;

% calcolo i polinomi di Lagrange associati ai nodi di Chebyshev
% i polinomi di Lagrange sono definiti come:
phi_1 = @(x) (x - x_v(2)) ./ (x_v(1) - x_v(2)) .* (x - x_v(3)) ./ (x_v(1) - x_v(3)) .* (x - x_v(4)) ./ (x_v(1) - x_v(4));
phi_2 = @(x) (x - x_v(1)) ./ (x_v(2) - x_v(1)) .* (x - x_v(3)) ./ (x_v(2) - x_v(3)) .* (x - x_v(4)) ./ (x_v(2) - x_v(4));
phi_3 = @(x) (x - x_v(1)) ./ (x_v(3) - x_v(1)) .* (x - x_v(2)) ./ (x_v(3) - x_v(2)) .* (x - x_v(4)) ./ (x_v(3) - x_v(4));
phi_4 = @(x) (x - x_v(1)) ./ (x_v(4) - x_v(1)) .* (x - x_v(2)) ./ (x_v(4) - x_v(2)) .* (x - x_v(3)) ./ (x_v(4) - x_v(3));

% calcolo la costante di lebesgue associata alla base di polinomi di lagrange

x_q = linspace(a, b, 1000); % punti di valutazione

leb = max(abs(phi_1(x_q)) + abs(phi_2(x_q)) + abs(phi_3(x_q)) + abs(phi_4(x_q)))
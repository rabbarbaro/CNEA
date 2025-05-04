clear
clc
close all

t_i = [0, 2, 5, 7, 8, 10]';
H_i = [400, 399.1, 397.5, 396.9, 396.4, 395.7]';

H0 = H_i(1);

% definisco il vettore di funzioni f (sui vari dati)
H = @(t, w) H0 - w(1) * (1 - exp(-w(2) * t));
% definisco il vettore dei residui (sui vari dati)
r = @(t, H_i, w) H_i - H(t, w);
% definisco lo jacobiano (derivo i residui, è una matrice)
J = @(t, H_i, w) [(1 - exp(-w(2) * t)), ...
                  w(1) * t .* (exp(-w(2) * t))];

w0 = [10, 0.1]';
toll = 1e-6;
maxit = 5;

% chiamiamo la funzione
[w_vect, it] = gauss_newton(t_i, H_i, r, J, w0, toll, maxit);
w_vect(:, end)

% plottiamo
plot(t_i, H_i, 'o', t_i, H(t_i, w_vect(:, end)));

t = linspace(0, 100);
plot(t_i, H_i, 'o', t, H(t, w_vect(:, end)));

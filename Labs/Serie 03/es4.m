clear
clc
close all

%% 1

% creiamo la funzione che implemennta l'algoritmo di iterazione QR
% vedi file qrbasic.m

%% 2

% definisco tolleranza come da testo
tol = 1e-10;
nmax = 1e3;

% costruiamo A1 con alpha = 30
A1 = [30 2 3 13
    5 11 10 8
    9 7 6 12
    4 14 15 1];
D1 = qrbasic(A1, tol, nmax);
fprintf("autovalori di A1:\n\t%f \n\t%f \n\t%f \n\t%f \n", D1)


% per verifica
% eig(A1)

% costruiamo A1 con alpha = 30
A2 = [-30 2 3 13
    5 11 10 8
    9 7 6 12
    4 14 15 1];
D2 = qrbasic(A2, tol, nmax);
fprintf("autovalori di A2:\n\t%f \n\t%f \n\t%f \n\t%f \n", D2)

% per verifica
% eig(A2)

% osseriviamo che la convergenza del metodo QR è più lenta per A2, dato che
% dipende dal massimo modulo dei rapporti degli autovalori

% come verifica:
v1 = abs(D1(2:end) ./ D1(1:end-1));
fprintf("massimo modulo dei rapporti degli autovalori di A1: %f\n", max(v1))
v2 = abs(D2(2:end) ./ D2(1:end-1));
fprintf("massimo modulo dei rapporti degli autovalori di A2: %f\n", max(v2))
clear
clc
close all

% inzializzo una matrice di Hilbert e una casuale
A = hilb(1000);
B = rand(1000);

% costruiamo i termini noti dei sistemi lineari in modo tale che le
% soluzioni siano vettori di uni
x_ex = ones(1000,1);
y_ex = ones(1000,1);
b = A * x_ex;
c = B * y_ex;

% risolviamo i sistemi con mldivide
x = A \ b;
y = B \ c;

% verifichiamo i numeri di condizionamento
cond(A)
cond(B)

% calcoliamo gli errori relativi
norm(x - x_ex) / norm(x_ex)
norm(y - y_ex) / norm(y_ex)

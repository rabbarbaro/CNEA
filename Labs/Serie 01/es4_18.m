clear
clc
close all

k = 6;

% dato che 5*k + k = 6*k è sempre pari
% creo due vettori che contengono i numeri pari e dispari nell'intervallo
x = (5*k + k):2:1000;
y = (5*k + k)+1:2:1000;

Sp = sum(x);
Sd = sum(y);
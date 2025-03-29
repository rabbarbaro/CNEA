clear
clc
close all

% definisco le funzioni
f1 = @(x) (x-1)^7;
f2 = @(x) x^7 - 7*x^6 + 21*x^5 - 35*x^4 + 35*x^3 - 21*x^2 + 7*x - 1;

x = 1.01;

% calcolo l'errore relativo percentuale
err_rel = 100 * abs(f1(x) - f2(x)) / abs(f1(x));
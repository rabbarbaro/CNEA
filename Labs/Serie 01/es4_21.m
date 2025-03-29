clear
clc
close all

% definisco k e l'intervallo di interi
k = 5;
x = 0:1000;

% verifico la non divisibilità e salvo in un vettore di valori logici
idx = mod(x, k+2) ~= 0;
% stampo la somma dei non divisibili
s = sum(x(idx));
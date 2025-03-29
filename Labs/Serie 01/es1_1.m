clear
clc
close all

% creo un vettore di potenze di 2, partendo da 2^0 fino a 2^10, passo 1
v1 = 2.^(0:10);

% creo un vettore (trasposto) di pi / intero da 1 a 10, passo 1
v2 = cos(pi./(1:10))';

% creo un vettore in cui divido 0.1 per potenze di 2 da 0 a 5, passo 1
v3 = 0.1 ./ (2.^(0:5));

% creo un vettore di zeri e riempio solo i componenti che mi interessano
% ci metto dentro una potenza di e da 1 a 10, passo 1
% ci sommo 6 + un multiplo di 5 da 5 a 45, cambio segno ogni volta
v4 = zeros (1,19);          % non sono sicuro sia necessario
v4 (1:2:19) = exp(1:10) - (-1).^(1:10) .* (6 + 5 .* (0:9));
clear
clc
close all

n = 1e6;

% di fatto è l'esercizio 4_4 dalla serie 01

% riempio una matrice con n coppie di numeri casuali fra 0 e 1
pnt = rand(2, n);
% calcolo il quadrato della distanza dal centro dei punti
pnt_dist = pnt(1, :).^2 + pnt (2, :).^2;
% verifico quali coppie sono nel cerchio e le conto
idx = pnt_dist <= 1;
pi_est = 4 * sum(idx) / n;
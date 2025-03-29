clear
clc
close all

n = 10;
N = 200;

% non sono riuscito a farlo in 1 o 2 istruzioni
% creo la matrice con la diagonale di 1, poi setto a mano gli altri 1
B = diag(ones(n, 1));
B(1:9:10, 2:9) = 1;
B(2:9, 1:9:10) = 1;

% la soluzione nota che si può ottere la stessa cosa da due prodotti
% vettore colonna per riga
% Bsol = diag(ones(1,10)) + [0,ones(1,8),0]' * [1,zeros(1,8),1] ...
%     + [1,zeros(1,8),1]' * [0,ones(1,8),0];
% isequal(B,Bsol);

C = diag(1:N) + diag(ones(N-1, 1), 1) + diag(ones(N-1, 1), -1) ...
    + diag(0.5 * ones(N-2, 1), 2) + 0.5 * diag(ones(N-2, 1), -2);

% per sapere se è corretta copio la soluzione e poi verifico
% Csol = diag([1:200]) + diag(ones(1,199),1) + diag(ones(1,199),-1) ...
%     + diag(0.5*ones(1,198),2) + diag(0.5*ones(1,198),-2);
% isequal(C,Csol);

% creo la matrice sommando le varie diagonali
% la prima sopradiagonale è un vettore di potenze di 3 da 0 a 8, passo 1
D = diag(20:-2:2) + diag(3.^(0:8), 1) + diag(0.5 * ones(8, 1), -2);
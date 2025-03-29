clear
clc
close all

n = 5;

% creo una matrice come richiesta sommando (sopra/sotto)diagonali
A = diag(2 .* (ones(n,1))) + diag(5 .* (ones(n-1,1)), -1) ...
    + diag(10 .* (ones(n-2,1)), 2) + diag(10 .* (ones(n-2,1)), -2) ...
    + diag(40 .* (ones(n-4,1)), 4) + diag(40 .* (ones(n-4,1)), -4);

% ottengo le sottomatrici richieste
A1 = A(1:3, 1:3);
A2 = A(1:2:5, [1 2 4]);
A3 = A(2:4, [1 3 4]);
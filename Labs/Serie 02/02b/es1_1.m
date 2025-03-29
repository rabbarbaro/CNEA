clear
clc
close all

%% 1

n = 100;
R_1 = 1;
R_2 = 2;

% costruiamo A
A = - diag(R_2 * ones(n, 1)) + diag(R_1 * ones(n-1, 1), -1);
A(1, :) = 1;

% stampo il numero di elementi non nulli e visualizzo la matrice
fprintf("Elementi non nulli di A: %d\n", nnz(A))
spy(A)

% rendo A sparsa
As = sparse(A);
% confronto la memoria occupata
fprintf("Dimensione di A piena: %d Bytes\n", getfield(whos('A'),'bytes'))
fprintf("Dimensione di A sparsa: %d Bytes\n", getfield(whos('As'),'bytes'))

%% 2

% calcolo la fattorizzazione LU di A
[L, U, P] = lu(A);

% visualizzo i pattern di sparsità di A, L, U
tiledlayout(2, 2);
nexttile
spy(A)
title("A")
nexttile
spy(L)
title("L")
nexttile
spy(U)
title("U")
nexttile
spy(P)
title("P")

% è avvenuto il fill-in

%% 3

% calcoliamo i raggi spettrali

% la matrice di iterazione g-s è B = I - T^-1 * A
P_gs = tril(A);
B_gs = eye(n) - P_gs \ A;       % inv(P_gs) * A
rho_B_gs = max(abs(eig(B_gs)));

% la matrice di iterazione j è B = I - D^-1 * A
P_j_inv = diag(1 ./ diag(A));
B_j = eye(n) - P_j_inv * A;
rho_B_j = max(abs(eig(B_j)));

fprintf("Raggio spettrale B G-S: %f\n", rho_B_gs)
fprintf("Raggio spettrale B J: %f\n", rho_B_j)

% il raggio spettrale di G-S è 1, il metodo non converge
% il raggio spettrale di J è < 1, il metodo converge

%% 4

% creiamo la funzione che applica il metodo di Jacobi
% vedi file jacobi.m

%% 5

% definisco il vettore dei termini noti
b = ones(n, 1);
b(1) = 2;
% definiamo i dati
x0 = zeros(n, 1);
toll = 1e-6;
nmax = 1e3;

% risolviamo con j
[x_j, k_j] = jacobi(A, b, x0, toll, nmax);

% verifichiamo che effettivamente sia la soluzione corretta
% x = A\b
% norm(x-x_j)/norm(x)
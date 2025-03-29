clear
clc
close all

% costruiamo la matrice (assicurata SDP) come da testo
n = 4;
A = 2 * diag(ones(n,1)) - diag(ones(n-1, 1), -1) - diag(ones(n-1, 1), 1);

%% 1

% creiamo la funzione che stima il numero di condizionamento della matrice
% usando il metodo delle potenze e il metodo delle potenze inverse per una
% matrice simmetrica e definita positiva
% per matrici SPD K(A) = K_2(A), ovvero se A SDP:
% sqrt(lambda_max(A'A)/lambda_min(A'A)) = lambda_max(A)/lambda_min(A)

% vedi file sdpcond.m

%% 2

% numero di iterazioni e tolleranza come da testo
tol = 1e-8;
nmax = 200;

% chiamo la funzione
K = sdpcond_sara(A, tol, nmax);

% per verifica
cond(A);

%% 3

% stima del numero di operazioni nella funzione
% O(2/3 * n^3)          LU fuori dal ciclo
% +
% O(iter * 8 * n^2)
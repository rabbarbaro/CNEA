clear
clc
close all

%% 1

% sapendo che gli autovalori di A: lambda_j(A) = 2 + 2*(cos(pi*j/(n+1)))
% per j = 1:n, lambda_1(A) > ··· > lambda_n(A)
% sapendo che il numero di condizionamento spettrale è dato dal rapporto
% fra l'autovalore di modulo massimo e l'autovalore di modulo minimo
% sapendo (dall'esercizio precedente) che A è SDN, quindi 
% la stima del numero di condizionamento spettrale K(A) di A in funzione di
% n per n-->+inf è:
% - j = 1: lim_n-->+inf lambda_1(A) = 2 + 2*(cos(pi/(n+1)))
%          lambda_1(A) = 2 + 2*(cos(0))
%          lambda_1(A) = 4
% - j = n: lim_n-->+inf lambda_n(A) = 2 + 2*(cos(pi*n/(n+1)))
%          lambda_n(A) = 2 + 2*(cos(pi))
%          lambda_n(A) = 0
%          NON POSSO DIVIDERE PER 0 --> espando in serie di Taylor il cos
%          cos(pi*x) = 1 − (pi*x)^2 / 2 + o(x^3)
%          cos(pi*n/(n+1)) = cos(pi*(1 - 1/(n+1)) = cos(pi - pi/(n+1))
%          cos(pi-x) = -cos(x)
%          cos(pi*n/(n+1)) = -cos(pi/(n+1)) = -(1 - (pi/(n+1))^2 / 2)
%          lim_n-->+inf lambda_n(A) = 2 + 2*(-1 + (pi/(n+1))^2 / 2)
%          lim_n-->+inf lambda_n(A) = (pi/(n+1))^2 = pi^2 / (n+1)^2
%          lambda_n(A) = pi^2 / n^2
% stima di K(A) = 4n^2 / pi^2

%% 2

% definisco il passo/numero di iterazioni N e inizializzo il vettore dei
% numeri di condizionamento spettrali
N = 10;
K = zeros(N, 1);

% per ogni iterazione
for ii = 1:N
    % definisco la dimensione del sistemo e genero A
    n = 10*ii;
    A = diag(2*(ones(n, 1))) ...
        - diag(ones(n-1, 1), -1)...
        - diag(ones(n-1, 1), 1);
    % calcolo il valore assoluto degli autovalori e divido il massimo per
    % il minimo (li salvo in un vettore per non doverli calcolare 2 volte)
    absEigA = abs(eig(A));
    K(ii) = max(absEigA)/min(absEigA);
end

% calcolo per ogni dimensione il numero di condizionamento spettrale
% stimato
n = 10 * (1:N);
K_est = 4 * n.^2 / pi^2;

% plotto
plot(n, K_est, n, K)
legend("K(A) stimato", "K(A)")

%% 3

% essendo la matrice tridiagonale (e lo sappiamo a priori) conviene
% utilizzare il metodo di Thomas O(8n), mentre la fattorizzazione LU è
% O(3/2 n)

%% 4

% per il calcolo del determinante posso ricordare che esso è una funzione
% moltiplicativa: det(A) = det(L)*det(U)
% sapendo a priori che la matrice è tridiagonale posso applicare la
% fattorizzazione di Thomas
% dato che sulla diagonale di L ci sono solo 1 il suo determinante è 1
% il determinante di A è quindi il determinante di U
% sulla diagonale di U ci sono gli elementi alpha, per calcolarli ho
% bisogno degli elementi delta
% devo quindi solo fattorizzare A e poi fare la produttoria degli elementi
% sulla diagonale di U
% sono:
% - (n-1) divisioni per calcolare delta
% - (n-1) sottrazioni e (n-1) moltiplicazioni per calcolare alpha
% - (n-1) moltiplicazioni per la produttoria
% il calcolo del determinante per matrici tridiagonali è O(4n), con 4n-4
% operazioni

% definisco una matrice tridiagonale arbitratia con soluzione esatta nota
% genero un vettore di termini noti
n = 10;
A = 2 * eye(n) - diag(ones(n-1, 1), 1) - diag(ones(n-1, 1), -1);
x_ex = ones(n, 1);
b = A * x_ex;

% calcolo con la funzione appena descritta
detA = detTriadiag(A);

% verifico che sia uguale al determinante calcolato da matlab
isequal(detTriadiag(A), det(A));


%% 5

% definisco matrice, vettore soluzione e termine noto come da testo
n = 100;
A = 2 * eye(n) - diag(ones(n-1, 1), 1) - diag(ones(n-1, 1), -1);
x_ex = ones(n, 1);
b = A * x_ex;

% costruisco il vettore c dato dal testo
c = rand(size(b));
c = c / norm(c);

% genero le perturbazioni su b
delta_b = 1e-6 * c;

% calcolo la soluzione perturbata e l'errore assoluto
x_pert = A \ (b + delta_b);
delta_x = x_pert - x_ex;

% calcolo l'errore relativo
err_rel = norm(delta_x) / norm(x_ex);

% essendo le perturbazioni su A nulle l'errore stimato è maggiorato come:
% err_rel_stim = cond(A) * norm(delta_b) / norm(b)

% calcolandolo verifichiamo che è una stima molto generosa
err_rel_stim = cond(A) * norm(delta_b) / norm(b);
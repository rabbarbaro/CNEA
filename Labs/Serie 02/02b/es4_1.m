clear
clc
close all

% genero la matrice A come dal testo
n = 300;
A = diag(6 * ones(n, 1)) - diag(2 * ones(n-1, 1), 1) ...
    - diag(2 * ones(n-1, 1), -1) + diag(ones(n-2, 1), 2) ...
    + diag(ones(n-2, 1), -2);

% data la soluzione esaatta estraggo il vettore dei termini noti
x_ex = ones(n, 1);
b = A * x_ex;

%% 1

% ricavo la matrice d'iterazione per il metodo di jacobi e calcolo il suo
% raggio spettrale
B_j = eye(n) - diag(1 ./ diag(A)) * A;
rho_B_j = max(abs(eig(B_j)));
fprintf("Raggio spettrale B_j: %f\n", rho_B_j)

% ricavo la matrice d'iterazione per il metodo di gauss-seidel e calcolo il
% suo raggio spettrale
P_gs = tril(A);
B_gs = eye(n) - P_gs \ A;
rho_B_gs = max(abs(eig(B_gs)));
fprintf("Raggio spettrale B_gs: %f\n", rho_B_gs)

% convergono entrambi, ma jacobi lo farà ad un rate molto molto lento
% come verifica:
% [x, k] = jacobi(A, b, zeros(n, 1), 1e-5, 1e5)
% per trovare una soluzione con toll = 1e-5 sono necessari ~7500 iterazioni

%% 2

% risolviamo con gauss-seidel, il raggio spettrale è molto inferiore
% rispetto a quello con jacobi

%% 3

% definisco la soluzione iniziale, la tolleranza e il numero massimo di
% iterazioni come dal testo
x0 = b;
toll = 1e-6;
nmax = 1000;

% risolvo con il metodo di gauss-seidel
[x_gs, k_gs] = gs(A, b, x0, toll, nmax);
fprintf("Iterazioni G-S: %d\n", k_gs)

% dalla teoria sappiamo che:
% res_norm = || b - A*x_k || / || b ||
% err_rel = || x_k - x_ex || / || x_ex ||
% err_rel_stima <= K(A) * res_norm

% calcolo l'errore relativo e il residuo normalizzato
err_rel = norm(x_ex - x_gs) / norm(x_ex);
res_norm = norm(b - A*x_gs) / norm(b);
fprintf("Errore relativo: %f\n", err_rel)
fprintf("Residuo normalizzato: %f\n", res_norm)

% dal residuo normalizzato e dal numero di condizionamento es
err_rel_stima = cond(A) * res_norm;
fprintf("Stima dell'errore relativo: <= %f\n", err_rel_stima)

%% 4

% dalla teoria sappiamo che:
% || e_k ||_A <= d^k || e_0 ||

% calcolo il numero di condizionamento di A, sapendo che essendo il metodo
% del gradente NON precondizionato P = I
d = (cond(A) - 1) / (cond(A) + 1);
% alternativamente: condA = max(eig(A)) / min(eig(A));
k = 20;

% calcolo l'errore in norma A alla prima iterazione (definizione) sapendo
% la soluzione esatta
err0_normA = sqrt((x_ex - x0)' * A * (x_ex - x0));

% calcolo la stima dell'errore in norma A dopo 20 iterazioni
err_normA_stima = d^k * err0_normA;
fprintf("Stima dell'errore in norma A per metodo del gradiente: <= %f\n", ...
    err_normA_stima)

%% 5

% creo il vettore dei valori beta
betav = 2:5;
% inizializzo il vettore con i d per ogni beta
dv = [];

% scorrendo tutti i beta
for beta = betav
    % genero la matrice precondizionatore P come da testo
    P = diag(beta * ones(n, 1)) - diag(ones(n-1, 1), 1) ...
    - diag(ones(n-1, 1), -1);
    % calcolo il numero di condizionamento di P^-1 * A
    condinvPA = cond(P\A);
    % alternativamente: condinvPA = max(eig(P\A)) / min(eig(P\A));
    % calcolo il valore d e lo "appendo" al vettore
    d = (condinvPA - 1 ) / (condinvPA + 1);
    dv = [dv, d];
end

% plotto per visualizzare meglio i risultati
plot(betav, dv)
% il metodo converge più rapidamente per beta = 4, ha il valore d minore
% dato che:
% || e_k ||_A <= d^k || e_0 ||
% essendo d < 0 l'errore tende a 0 per k che tende all'infinito, e ci
% tenderà tanto più velocemente tanto più è piccolo d

%% 6

% sappiamo che sicuramente (in aritmetica esatta) converge alla soluzione
% esatta al più in n = 300 iterazioni, per la tolleranza eichiesta 100 è
% più che sufficiente
nmax_cg = 100;

% chiamo la funzione gradiente coniugato precondizionato in cui la matrice
% di precondizionamento è l'identità (non precondizionato)
[x_cg, ~, ~, k_pcg] = pcg(A, b, toll, nmax_cg, eye(n), eye(n), x0);
fprintf("Iterazioni GC: %d\n", k_pcg)

% conoscendo la soluzione esatta calcolo la norma A dell'errore
% direttamente
err_A_pcg = sqrt((x_cg - x_ex)' * A * (x_cg - x_ex));
fprintf("Errore in norma A per metodo del gradiente coniugato: %f\n", ...
    err_A_pcg)
clear
clc
close all

%% 1

% generiamo la matrice, il termine noto e la soluzione iniziale
n = 50;
A = diag(4 * ones(n, 1)) - diag(ones(n-1, 1), 1) ...
    - diag(ones(n-1, 1), -1) - diag(ones(n-2, 1), 2) ...
    - diag(ones(n-2, 1), -2);
b = 0.2 * ones(n,1);
x0 = zeros(n,1);

toll = 1e-5;

%% 2

if A == A'
    if eig(A) > 0
        disp("A SDP")
    end
end

condA = max(eig(A)) / min(eig(A));
fprintf("Numero di condizionamento di A: %f\n", condA)
% verifica
% cond(A)

%% 3

% creiamo la funzione che applica il metodo di Richardson
% vedi file richardson.m

%% 4

% per avere richardson non precondizionato pongo P = I
P = eye(n);
nmax = 1e4;

% metodo statico, assegno alpha (casi 1, 2, opt)
% calcolo la matrice di iterazione per ogni valore di alpha e il suo raggio
% spettrale, poi calcolo la soluzione con quell'alpha

alpha_1 = 0.2;
B_alpha_1 = eye(n) - alpha_1 * A; 
rho_B_alpha_1 = max(abs(eig(B_alpha_1)));

[x_1, k_1] = richardson(A, b, P, x0, toll, nmax, alpha_1);

alpha_2 = 0.33;
B_alpha_2 = eye(n) - alpha_2 * A; 
rho_B_alpha_2 = max(abs(eig(B_alpha_2)));

[x_2, k_2] = richardson(A, b, P, x0, toll, nmax, alpha_2);

% dato che P = I, P^-1 * A = A
lambda_max = max(eig(A));
lambda_min = min(eig(A));
alpha_opt = 2 / (lambda_min + lambda_max);
B_alpha_opt = eye(n) - alpha_opt * A;
rho_B_alpha_opt = max(abs(eig(B_alpha_opt)));

[x_opt, k_opt] = richardson(A, b, P, x0, toll, nmax, alpha_opt);

fprintf("Iterazioni con alpha_1 = 0.2: %d\n", k_1)
fprintf("Iterazioni con alpha_2 = 0.33: %d\n", k_2)
fprintf("Iterazioni con alpha_opt: %d\n", k_opt)
fprintf("Raggio spettrale con alpha_1 = 0.2: %f\n", rho_B_alpha_1)
fprintf("Raggio spettrale con alpha_2 = 0.33: %f\n", rho_B_alpha_2)
fprintf("Raggio spettrale con alpha_opt: %f\n", rho_B_alpha_opt)

% usando alpha_2 = 0.33 NON converge, il raggio spettrale è > 1
% usando alpha_opt avremo meno iterazioni rispetto a usare alpha_2

% per curiosità vediamo come si comporta il metodo dinamico (gradiente)
[x_dyn, k_dyn] = richardson(A, b, P, x0, toll, nmax);
fprintf("Iterazioni metodo dinamico (gradiente): %d\n", k_dyn)
% in questo caso converge leggermente più lentamente del metodo statico con
% alpha_opt

%% 5

% definisco il precondizionatore e alpha come de testo
P = tril(A);
alpha = 1;

% calcolo il raggio spettrale di B
B_r = eye(n) - alpha * (P \ A);
rho_B_r = max(abs(eig(B_r)));

% calcolo con richardson
[x_r, k_r] = richardson(A, b, P, x0, toll, nmax, alpha);

% noto che, essendo P = tril(A) e alpha = 1 costante, in questo caso il
% metodo di richardson coincide esattamente con il metodo di g-s
[x_gs, k_gs] = gs(A, b, x0, toll, nmax);

fprintf("Iterazioni con Richarson con P = tril(A) e alpha = 1: %d\n", k_r)
fprintf("Iterazioni con Gauss-Seidel: %d\n", k_gs)
fprintf("Raggio spettrale con P = tril(A) e alpha = 1: %f\n", rho_B_r)

%% 6

% costruisco il precondizionatore come nel testo
P = diag(2 * ones(n, 1)) - diag(ones(n-1, 1), 1) - diag(ones(n-1, 1), -1);

% verifico che P sia SDP
if P == P'
    if eig(P) > 0
        disp("P SPD")
    end
end

% calcolo il numero di condizionamento di P^-1 * A
condinvPA = cond(P \ A);
fprintf("Numero di condizionamento di P^-1 * A: %f\n", condinvPA)

% calcolo alpha ottimale e il raggio spettrale della matrice B
lambda_max_invPA = max(eig(P\A));
lambda_min_invPA = min(eig(P\A));
alpha_r_opt = 2 / (lambda_min_invPA + lambda_max_invPA);
B_r_opt = eye(n) - alpha_r_opt * (P \ A);
rho_B_r_opt = max(abs(eig(B_r_opt)));

% risolvo con richardson stazionario
[x_r_opt, k_r_opt] = richardson(A, b, P, x0, toll, nmax, alpha_r_opt);

% risolvo con il metodo del gradiente precondizionato
[x_grad, k_grad] = richardson(A, b, P, x0, toll, nmax);

fprintf("Raggio spettrale B: %f\n", rho_B_r_opt)
fprintf("Iterazioni con Richardson con alpha_r_opt: %d\n", k_r_opt)
fprintf("Iterazioni con gradiente precondizionato: %d\n", k_grad)

% in questo caso il metodo del gradiente converge più velocemente del
% metodo stazionario
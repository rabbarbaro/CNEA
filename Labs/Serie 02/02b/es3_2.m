clear
clc
close all

%% 1

% definisco il vettore con le dimensioni dei sistemi
N = 5:20;

% definisco tolleranza e iterazioni massime come da testo
toll = 1e-6;
nmax = 5e3;

% inizializzo i vettori delle iterazioni per il metodo del gradiente non
% precondizionato, precondizionato, del numero di condizionamento di A, del
% numero di condizionamento di P^-1 * A
it_np = zeros(1, length(N));
it_p = zeros(1, length(N));
k_A = zeros(1, length(N));
k_PinvA = zeros(1, length(N));

for n = N
    % genero A, b e x0 come da testo
    A = diag(4 * ones(n, 1)) + diag(ones(n-1, 1), 1) ...
        + diag(ones(n-1, 1), -1) + diag(2*ones(n-2, 1), 2) ...
        + diag(2*ones(n-2, 1), -2);
    b = ones(n, 1);
    x0 = zeros(n,1);
    
    % calcolo con il metodo del gradiente non precondizionato
    [x_np, k_np] = richardson(A, b, eye(n), x0, toll, nmax);

    % % calcolo con il metodo del gradiente precondizionato
    P = tril(A);
    [x_p, k_p] = richardson(A, b, P, x0, toll, nmax);

    % salvo per questa dimensione il numero delle iterazioni
    it_np(n - 4) = k_np;
    it_p(n - 4) = k_p;
    % salvo anche i numeri di condizionamento
    k_A(n - 4) = cond(A);
    k_PinvA(n - 4) = cond(P\A);
end

%% 2 + 3

% rappresento in scala logaritmica il numero di iterazioni 
tiledlayout(1, 2)
nexttile
semilogy(N, it_np, N, it_p)
grid on
title("Iterazioni G VS G precondizionato")
legend("Iterazioni gradiente non precondizionato", ...
    "Iterazioni gradiente precondizionato")
nexttile
plot(N, k_A, N, k_PinvA)
grid on
title("Numero di condizionamento A e P^-1 * A")
legend("Numero condizionamento A", ...
    "Numero condizionamento P^-1*A")

%% 4

% inizializzo il vettore delle iterazioni
it_cg = zeros(1, length(N));

for n = N
    % definisco A, b, x0 per ogni dimensione
    A = diag(4 * ones(n, 1)) + diag(ones(n-1, 1), 1) ...
        + diag(ones(n-1, 1), -1) + diag(2*ones(n-2, 1), 2) ...
        + diag(2*ones(n-2, 1), -2);
    b = ones(n, 1);
    x0 = zeros(n,1);
    
    % chiamo la funzione per il gradiente coniugato precondizionato usando
    % come precondizionatore l'identità (non precondizionato)
    [x_cg, ~, ~, k_cg] = pcg(A, b, toll, nmax, eye(n), eye(n), x0);

    % salvo
    it_cg(n - 4) = k_cg;
end

figure
semilogy(N, it_np, N, it_p, N, it_cg)
grid on
legend("Iterazioni gradiente non precondizionato", ...
    "Iterazioni gradiente precondizionato", ...
    "Iterazioni gradiente coniugato non precondizionato")
title("Confronto numero interazioni fra G, GP, CG")
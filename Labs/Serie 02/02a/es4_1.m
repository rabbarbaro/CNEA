clear
clc
close all

% prealloco la matrice contenente i tempi di esecuzione
T = zeros(4, 10);

% definisco il numero di iterazioni
iter = 10;

% faccio da 1 a 10, in modo da avere subito ii (dovrò poi moltiplicare per
% 500 per ottenere il passo richiesto su n
for ii = 1:iter
    % inizializzo i dati
    K = 100;
    L = 20;
    n = 500 * ii;
    
    % genero la matrice A
    A = diag(-2*(ones(n, 1))) + ...
        diag(ones(n-1, 1), -1) + ...
        diag(ones(n-1, 1), 1);
    A = K*A;
    
    % definisco il vettore f
    f = zeros(n,1);
    f(end) = -K * L;
    
    % risolvo con fattorizzazione LU e poi i due sistemi lineari
    % triangolari back-to-back con mldivide
    tic
    [LL, UU, PP] = lu(A);
    x_lu = UU \ (LL \ (PP * f));
    T(1,ii) = toc;
    
    % risolvo con Thomas (sapendo la tridiagonalità a priori)
    tic
    [LL, UU, x] = thomas(A, f);
    T(2,ii) = toc;

    % risolvo con Thomas (verificando la tridiagonalità)
    tic
    [LL, UU, x] = thomas(A, f, 1);
    T(3,ii) = toc;
    
    % risolvo con Cholesky
    % dato che la matrice A è SDN, e Cholesky si applica solo a matrici
    % SDP, calcolando l'opposto sia di A che di f otterremo il risultato
    % corretto e potremo utilizzare questo algoritmo
    % per "correttezza" eseguiamo questa operazione fuori dal tic-toc
    A_chol = -A;
    f_chol = -f;
    % calcolo la fattorizzazione di Cholesky con il comando di matlab e
    % risolvo i due sistemi triangolari back-to-back con mldivide
    tic
    H = chol(A_chol);
    x_chol = H \ (H' \ f_chol);               % -A = h' * H
    T(4,ii) = toc;
end

%% plot dei tempi

dim = 500 * [1:iter];
plot(dim, T(1, :), '-o')
hold on
plot(dim, T(2, :), '-o')
plot(dim, T(3, :), '-o')
plot(dim, T(4, :), '-o')
legend('LU', 'Thomas', 'Thomas w/ verifica', 'Cholesky')
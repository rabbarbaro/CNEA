function K = sdpcond_mod(A, tol, nmax)

% K = sdpcond(A, tol, nmax)
% Metodo per il calcolo del numero di condizionamento spettrale di matrici
% simmetriche e definite positive
% IN
%   - A: matrice sistema
%   - tol: tolleranza criterio d'arresto (differenza iterate successive)
%   - nmax: numero massimo di iterazioni
% OUT
%   - K: numero di condizionamento spettrale di A

%% inizializzazione

% calcoliamo un'unica volta la fattorizzazione LU
[L, U, P] = lu(A);                                                          % O(2/3 n^3) (operazione più costosa fuori dal ciclo)

% generiamo un vettore di partenza casuale e inizializziamo le iterazioni a
% zero
n = size(A, 1);
x0 = ones(n, 1);
iter = 0;

% normalizzo il vettore di partenza per il metodo diretto
y = x0 / norm(x0);
% il vettore per il metodo inverso è lo stesso, non ricalcolo x0 / norm(x0)
y_mu = y;

% calcolo la stima iniziale dell'autovalore massimo
lambda = y' * A * y;
% la stima inizale dell'autovettore minimo è la stessa, non ricalcolo
% y' * A * y;
mu = lambda;

% inizializziamo il numero di condizionamento a 1 (lambda = mu)
K0 = 1;

% poniamo l'errore iniziale arbitrario maggiore della tolleranza
err = tol + 1;

% usiamo l'incremento tra iterate successive come stimatore dell'errore
while err > tol && iter < nmax
    iter = iter + 1;

    % calcoliamo la stima successiva dell'autovalore massimo
    x = A * y;                                                              % 2n^2-n ops
    y = x / norm(x);                                                        % 2n+1 ops 
    lambda = y' * A * y;                                                    % 2n^2+n-1 ops

    % calcoliamo la stima successiva dell'autovalore minimo
    z = fwsub(L, P*y_mu);                                                   % n^2 ops
    x = bksub(U, z);                                                        % n^2 ops
    y_mu = x / norm(x);                                                     % 2n+1 ops
    mu = y_mu' * A * y_mu;                                                  % 2n^2+n-1 ops

    % aggiorniamo il numero di condizionamento spettrale
    K = lambda/mu;

    % calcoliamo l'errore (non normalizzato)
    err = abs(K- K0);

    % salviamo il valore corrente della stima per aggiornare l'errore
    K0 = K;
end

end
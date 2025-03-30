function [lambda, x, iter] = eigpower(A, tol, nmax, x0)

% [lambda, x, iter] = eigpower(A, tol, nmax, x0)
% Metodo delle potenze per approssimare l'autovalore di modulo massimo e il
% suo autovettore corrispondente, Ax = lambdax
% IN
%   - A: matrice quadrata (n x n)
%   - tol: tolleranza sul criterio d'arresto (differenza iterate relativa)
%   - nmax: numero massimo di iterazioni
%   - x0: iterata iniziale per l'autovettore, vettore colonna
% OUT
%   - lambda: approssimazione dell'autovalore di modulo massimo
%   - x: approssimazione dell'autovettore corrispondente a lambda (non
%       normalizzato)
%   - iter: numero di iterazioni effettuate

%% verifica input

% dimensioni
[n,m] = size(A);
if n ~= m
    error('Dimensioni incompatibili (solo per matrici quadrate)');
end

%% inizializzazione

% poniamo le iterazioni a 0, normalizziamo l'iterata iniziale e calcoliamo
% la stima iniziale dell'autovalore con quoziente di rayleigh
iter = 0;
y = x0 / norm(x0);
lambda0 = y' * A * y;

% poniamo errore iniziale arbitrario maggiore della tolleranza
err = tol + 1;

%% algoritmo

% finché non rispettiamo la tolleranza, non superiamo le iterazioni massime
% e l'autovalore non è nullo
while err > tol && iter < nmax && abs(lambda0) ~= 0
    iter = iter + 1;
    % aggiorniamo l'autovettore e lo normalizziamo
    x = A * y;
    y = x / norm(x);
    % aggiorniamo autovalore con quoziente di rayleigh
    lambda = y' * A * y;
    % calcoliamo la stima dell'errore normalizzata su lambda0 (divideremmo
    % per lambda, ma dobbiamo verificare che sia diverso da 0)
    err = abs(lambda - lambda0) / abs(lambda0);
    % salviamo il valore corrente della stima per aggiornare l'errore
    lambda0 = lambda;
end

end
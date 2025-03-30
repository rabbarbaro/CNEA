function D = qrbasic(A, tol, nmax)
% D = qrbasic(A, tol, nmax)
% Calcola con il metodo delle iterazioni QR tutti gli autovalori della
% matrice A (convergenza non garantita in generale)
% IN
%   - A: matrice quadrata (n x n)
%   - tol: tolleranza sul criterio d'arresto (massimo modulo nella matrice
%       triangolare inferiore, diagonale esclusa)       
%   - nmax: numero massimo di iterazioni
% OUT
%   - D: approssimazione degli autovalori di A, vettore colonna 

%% verifica input

% dimensioni
[n, m] = size(A);
if n ~= m
    error('La matrice deve essere quadrata')
end

%% inizializzaione

% azzeriamo il contatore delle iterazioni
k = 0;

% abbiamo due scelte come inizio:
% poniamo uno stimatore dell'errore arbitrario maggiore della tolleranza,
% costa poco ma se la matrice è già diagonale facciamo un'iterazione in
% più inutile
est = tol + 1;
% calcoliamo lo stimatore, se è già diagonale usciamo subito, ma costa di
% più da calcolare
% est = max(max(abs(tril(A, -1))));

%% algoritmo

% finché non superiamo le iterazioni massime e lo stimatore è maggiore
% della tolleranza
while k < nmax && est > tol
    k = k + 1;
    % calcolo la fattorizzazione QR
    [Q, R] = qr(A);
    % aggiorno la matrice A
    A = R * Q;
    % calcolo lo stimatore dell'errore come il massimo modulo nella matrice
    % triangolare inferiore, diagonale esclusa
    est = max(max(abs(tril(A, -1))));
end

% estraggo il vettore con gli autovalori di A dalla diagonale
D = diag(A);

if k ~= nmax
    fprintf("qrbasic converge in %d iterazioni \n", k)
else
    fprintf("qrbasic NON converge in %d iterazioni \n", k)
end

end
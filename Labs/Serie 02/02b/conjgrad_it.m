function [xk, k] = conjgrad_it(A, b, x0, nmax, toll)

% [x, k] = conjgrad:it(A, b, P, x0, nmax, toll)
% Metodo del gradiente coniugato non precondizionato per risolvere
% iterativamente sistemi lineari Ax = b
% IN
% - A: matrice del sistema
% - b: vettore termine noto
% - x0: guess iniziale
% - nmax: nmax: massimo numero di iterazioni
% - toll: tolleranza sul residuo normalizzato
% OUT
% - x: soluzione ottenuta
% - k: numero di iterazioni effettuate

%% verifica input

n = length(b);
k = 0;

if ((size(A,1) ~= n) || (size(A,2) ~= n) || (length(x0) ~= n))
  error('Dimensioni incompatibili')
end

%% inizializzazione

% soluzione, residuo, direzione e residuo normalizzato all'iterata 0
x = x0;
r = b - A*x;
p = r;
err = norm(r) / norm(b);

xk = x0;

%% algoritmo

% finché il residuo normalizzato è superiore alla tolleranza e non abbiamo
% raggiunto il limite di iterazioni
% seguo l'algoritmo
while err > toll && k < nmax
    k = k + 1;
    alpha = (p' * r) / (p' * A * p);
    x = x + alpha * p;
    r = r - alpha * A * p;
    beta = (p' * A * r) / (p' * A * p);
    p = r - beta * p;

    err = norm(r) / norm(b);
    xk = [xk x];
end

end
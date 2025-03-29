function [x, k] = jacobi(A, b, x0, toll, nmax)

% [x, k] = jacobi(A, b, x0, toll, nmax)
% Metodo di Jacobi per risolvere iterativamente sistemi lineari Ax = b
% IN
% - A: matrice del sistema
% - b: vettore termine noto
% - x0: vettore iniziale
% - toll: tolleranza sul residuo normalizzato
% - nmax: massimo numero di iterazioni
% OUT
% - x: soluzione ottenuta
% - k: numero di iterazioni effettuate

%% verifica degli input

n = length(b);
k = 0;

% verifica delle dimensioni
if ((size(A,1) ~= n) || (size(A,2) ~= n) || (length(x0) ~= n))
  error('dimensioni incompatibili')
end

% controllo che gli elementi diagonali non siano nulli
if (prod(diag(A)) == 0)
  error('elementi diagonali nulli')
end

%% inizializzazione

% P = D = diag(A)
D_inv = diag(1 ./ diag(A));
% soluzione, residuo e residuo normalizzato all'iterata 0
x = x0;
r = b - A*x;
err = norm(r) / norm(b);

%% algoritmo

% finché il residuo normalizzato è superiore alla tolleranza e non abbiamo
% raggiunto il limite di iterazioni
while err > toll && k < nmax
    k = k + 1;
    % risolvo il sistema Pz = r
    z = D_inv * r;
    % aggiorno x ed r
    x = x + z;
    r = b - A*x;
    % aggiorno il residuo normalizzato
    err = norm(r) / norm(b);
end

end
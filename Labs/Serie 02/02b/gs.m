function [x, k] = gs(A, b, x0, toll, nmax)

% [x, k] = gs(A, b, x0, toll, nmax)
% Metodo di Gauss-Seidel per risolvere iterativamente sistemi lineari
% Ax = b
% IN
% - A: matrice del sistema
% - b: vettore termine noto
% - x0: vettore iniziale
% - toll: tolleranza sul residuo normalizzato
% - nmax: massimo numero di iterazioni
% OUT
% - x: soluzione ottenuta
% - k: numero di iterazioni effettuate

%% verifica input

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

T = tril(A);
% soluzione, residuo e residuo normalizzato all'iterata 0
x = x0;
r = b - A*x;
err = norm(r) / norm(b);

% finché il residuo normalizzato è superiore alla tolleranza e non abbiamo
% raggiunto il limite di iterazioni
while err > toll && k < nmax
    k = k + 1;
    % risolvo il sistema Tz = r (con fwsub perché so essere triangolare)
    z = fwsub(T, r);
    % aggiorno x ed r
    x = x + z;
    r = b - A*x;
    % aggiorno il residuo normalizzato
    err = norm(r) / norm(b);
end

end
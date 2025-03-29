function [x, k] = richardson(A, b, P, x0, toll, nmax, alpha)

% [x, k] = richardson(A, b, P, x0, toll, nmax, (opt) alpha)
% Metodo di Richardson stazionario precondizionato o dinamico
% precondizionato (gradiente precondizionato) per risolvere iterativamente
% sistemi lineari Ax = b
% IN
% - A: matrice del sistema
% - b: vettore termine noto
% - P: matrice precondizionatore
% - x0: guess iniziale
% - toll: tolleranza sul residuo normalizzato
% - nmax: nmax: massimo numero di iterazioni
% - (opt) alpha: opzionale, parametro di accelerazione
%       - se non assegnato si considera il metodo dinamico (gradiente
%         precondizionato)
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
    z = P \ r;
    % se NON mi è stato assegnato alpha calcolo a ogni iterazione alpha
    if nargin == 6
        alpha = (z' * r) / (z' * A * z);
    end
    % aggiorno x ed r
    x = x + alpha * z;
    r = b - A*x;
    % aggiorno il residuo normalizzato
    err = norm(r) / norm(b);
end

end
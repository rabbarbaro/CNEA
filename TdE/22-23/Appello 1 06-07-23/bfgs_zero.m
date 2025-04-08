function [xvect, it] = bfgs_zero(F, B0, x0, toll, nmax)

% [xvect, it] = bfgs_zero(F, B0, x0, toll, nmax)
% Algoritmo BFGS per l'approssimazione di zeri per sistemi di equazioni non
% lineari
% IN
%   - F: funzione vettoriale R^n --> R^n
%   - B0: approssimazione della jacobiana di F
%   - x0: guess iniziale
%   - toll: tolleranza
%   - nmax: numero massimo di iterazioni
% OUT
%   - xvect: matrice che su ogni colonna contiene la soluzione a ogni
%       iterazione (l'ultima colonna è la soluzione)
%   - it: numero di iterazioni

% inizializzo variabili per ciclo iterativo
it = 0;
err = toll+1;

% inizializzo la matrice soluzione
xvect = x0;

% definisco la prima approssimazione della jacobiana
B = B0;

% ciclo while per calcolo x^{it+1} (it = 0, 1, ...)
while it < nmax && err > toll
    
    % risolvo il sistema per trovare la differenza
    delta = B \ (-F(x0));

    % calcolo il nuovo valore di x
    x_new = x0 + delta;

    % calcolo y
    y = F(x_new) - F(x0);

    % aggiorno B
    B = B + (y * y') / (y' * delta) - (((B * delta) * (B * delta)')) / (delta' * (B * delta));

    % Aggiorno quantità per metodo iterativo
    it = it+1;
    err = norm(F(x_new));
    x0 = x_new;
    xvect = [xvect, x_new];
end

% printo info utili
if err <= toll
     fprintf("Il metodo BFGS converge in %d iterazioni alllo zero alpha\n", it);
else
     fprintf("Il metodo BFGS non converge in %d iterazioni. \n", it)
end

end
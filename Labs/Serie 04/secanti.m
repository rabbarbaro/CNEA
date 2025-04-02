function [xvect, it] = secanti(xStart, nmax, toll, fun)

% [xvect, it] = secanti(x0, nmax, toll, fun)
% Metodo delle secanti per l'approssimazione degli zeri della funzione fun.
% Criterio d'arresto basato sul controllo della differenza tra due iterate
% successive.
% IN
%   - xStart: punti di partenza
%   - nmax: numero massimo di iterazioni
%   - toll: tolleranza sul test d'arresto
%   - fun: funzione come anonymus/inline
% OUT
%   - xvect: vettore (riga) con tutte le iterate calcolate (l'ultima è
%       la soluzione)
%   - it: iterazioni effettuate

% per entrare nel ciclo, inizializziamo errore, iterazioni e soluzioni
err = toll + 1;
it = 1;
xvect = xStart;

% finché non abbiamo raggiunto il numero di iterazioni massimo e l'errore
% è meggiore della tolleranza
while it < nmax && err > toll
    % calcolo q
    q = (fun(xvect(end)) - fun(xvect(end - 1))) / (xvect(end) - xvect(end - 1));
    % aggiorno la soluzione
    xvect = [xvect, xvect(end) - fun(xvect(end))/q];
    % calcolo stimatore errore
    err = abs(xvect(end)- xvect(end - 1));
    it = it + 1;
end

if (it < nmax)
    fprintf("Convergenza al passo %d \n", it);
else
    fprintf("Raggiunto il numero massimo di passi %d \n", it);
end
fprintf("Radice calcolata: %-12.8f \n", xvect(end))
function [xmatrix, it] = newton_sistemi( x0, nmax, toll, F, J)

% Metodo di Newton per la ricerca della soluzione di sistemi di
% equazioni non lineari. 
% Criterio d'arresto basato sul controllo della norma del residuo.

% Parametri di ingresso:
%
%   x0: Punto di partenza (vettore)
%   nmax: Numero massimo di iterazioni
%   toll: Tolleranza sul test d'arresto
%   fun: sistema di equazioni lineari
%   jac: matrice jacobiana del sistema


% Parametri di uscita:
%
%   xmatrix: matrice contenete tutte le iterate (che sono dei vettori 2*1)
%            e l'ultima componente e' la soluzione
%   it: Iterazioni effettuate


err = toll + 1;
it = 0;
xmatrix = x0;
x = x0;

while it < nmax && err> toll
   
    sigma = - J( x ) \ F( x );
    xnew = x + sigma;
    
    err = norm( F(xnew) );
    
    xmatrix = [xmatrix, xnew];
    it = it + 1;
    x = xnew;
   
end

% controllo il motivo per cui sono uscita dal ciclo e stampo risultati
if (it < nmax)
    fprintf(' Convergenza al passo k : %d \n', it);
else
    fprintf(' E` stato raggiunto il numero massimo di passi k : %d \n', it);
end
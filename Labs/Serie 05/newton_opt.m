function [xvect, it] = newton_opt(GradPhi, HessPhi, x0, toll, nmax)

% [xvect, it] = newton_opt(GradPhi, HessPhi, x0, toll, nmax)
% Algoritmo di Newton per la risoluzione di problemi di ottimizzazione per
% funzionali in 2-dimensioni. Consideriamo il metodo di Newton esatto,
% quindi fissiamo \alpha_k = 1
% IN
%   - GradPhi: gradiente del funzionale
%   - HessPhi: hessiano del funzionale 
%   - x0: guess iniziale
%   - toll: tolleranza
%   - nmax: numero massimo di iterazioni
% OUT
%   - xvect: matrice che su ogni colonna contiene la soluzione a ogni
%       iterazione (l'ultima colonna è la soluzione)
%   - it: numero di iterazioni

% inizializzo variabili per ciclo iterativo
it = 0;
err = toll + 1;

% inizializzo la matrice soluzione
xvect = x0;

% dalcolo la direzione iniziale di discesa
d = - HessPhi(x0(1),x0(2)) \ GradPhi(x0(1),x0(2));

% definiamo lo step-size alpha_k
alpha_k = 1;

% ciclo while per calcolo x^{it+1} (it = 0, 1, ...)
while it < nmax && err > toll
    
    % calcolo il nuovo valore di x
    x_new = x0 + alpha_k*d;
    
    % aggiorno la direzione di discesa
    d = - HessPhi(x_new(1), x_new(2)) \ GradPhi(x_new(1), x_new(2));

    % aggiorno quantità per metodo iterativo
    it = it + 1;
    err = norm(GradPhi(x_new(1), x_new(2)));
    x0 = x_new;
    xvect = [xvect, x_new];
end

if err <= toll
     fprintf("Il metodo di Newton converge in %d iterazioni al punto di minimo x = (%f,%f)^T \n", it, xvect(1,end), xvect(2,end));
else
     fprintf("Il metodo di Newton non converge in %d iterazioni. \n", it)
end

end
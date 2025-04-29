function [w_vect, it] = gauss_newton(x, y, r_fun, J_fun, w0, toll, maxit)

% [w_vect, it] = gauss_newton(x, y, r_fun, J_fun, w0, toll, maxit)
% Algoritmo di Gauss-Newton per il calcolo del minimo di una funzione non
% lineare con metodo dei minimi quadrati non lineari.
% IN
%   - x: vettore colonna di dimensione n x 1 (n = numero di dati)
%   - y: vettore colonna di dimensione n x 1
%   - r_fun: anonymus function del residuo (n x 1)
%   - J_fun: anonymus function dello jacobiano (n x m)
%         (m = numero di parametri)
%   - w0: vettore colonna di dimensione m x 1 (guess iniziale)
%   - toll: tolleranza per il criterio di arresto (scelta dell'errore)
%   - maxit: numero massimo di iterazioni
% OUT
%   - w_vect: matrice colonna di dimensione m x (it + 1) con i valori
%           della soluzione ad ogni iterazione (e la guess iniziale)
%   - it: numero di iterazioni effettuate

% inizializzo variabili per ciclo iterativo
it = 0;

% inizializzo la matrice soluzione con la guess iniziale
w_vect = w0;

% calcolo jacobiano e residuo iniziali
J = J_fun(x, y, w0);
r = r_fun(x, y, w0);

% dalcolo la direzione iniziale di discesa
d = - (J' * J) \ (J' * r);

% uso la norma della direzione di discesa come errore per il criterio di
% arresto
err = norm(d);

while it < maxit && err > toll
    % calcolo il nuovo valore di w
    w_new = w0 + d;
    
    % aggiorno jacobiano e residuo
    J = J_fun(x, y, w_new);
    r = r_fun(x, y, w_new);
    
    % aggiorno la direzione di discesa
    d = - (J' * J) \ (J' * r);

    % aggiorno quantità per metodo iterativo
    it = it + 1;
    err = norm(d);
    w0 = w_new;
    w_vect = [w_vect, w_new];
end

if err <= toll
     fprintf("Il metodo di Gauss-Newton converge in %d iterazioni\n", it)
else
     fprintf("Il metodo di Gauss-Newton non converge in %d iterazioni. \n", it)
end

end
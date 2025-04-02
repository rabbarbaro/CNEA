function [xvect, it] = bfgs(Phi, GradPhi, x0, toll, nmax)

% [xvect, it] = bfgs(Phi, GradPhi, x0, toll, nmax)
% Algoritmo BFGS per la risoluzione di problemi di ottimizzazione per
% funzionali in 2-dimensioni. Consideriamo il metodo di backtracking per
% aggiornare lo step-size alpha_k
% IN
%   - Phi: funzionale (non necessario per il BFGS ma per il backtracking)
%   - GradPhi: gradiente del funzionale
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

% calcolo la direzione iniziale di discesa
I = eye(length(x0));
Bk = I;
d = - Bk * GradPhi(x0(1), x0(2));

% parametri backtracking
c1_bt = 1e-4;
rho_bt = 0.5;
nmax_bt = 10;

% ciclo while per calcolo x^{it+1} (it = 0, 1, ...)
while it < nmax && err > toll
    
    % calcolo il nuovo valore di alpha_k con backtracking
    [alpha_k, ~] = backtracking(Phi, GradPhi, x0, d, c1_bt, rho_bt, nmax_bt);    
    
    % calcolo il nuovo valore di x
    x_new = x0 + alpha_k*d;

    % calcolo Bk (approssimazione inversa Hessiano di Phi)
    delta_k = x_new - x0;
    s_k = GradPhi(x_new(1), x_new(2)) - GradPhi(x0(1), x0(2));  
    rho_k = 1 / (delta_k' * s_k);    
    Bk = (I - rho_k*delta_k*s_k')*Bk*(I - rho_k*s_k*delta_k') + (rho_k*delta_k)*delta_k';    
    
    % Aggiorno la direzione di discesa
    d = - Bk*GradPhi(x_new(1), x_new(2));

    % Aggiorno quantità per metodo iterativo
    it    = it+1;
    err   = norm(GradPhi(x_new(1), x_new(2)));
    x0    = x_new;
    xvect = [xvect, x_new];
    
end

if err <= toll
     fprintf("Il metodo BFGS converge in %d iterazioni al punto di minimo x = (%f,%f)^T \n", it, xvect(1,end), xvect(2,end));
else
     fprintf("Il metodo BFGS non converge in %d iterazioni. \n", it)
end

end
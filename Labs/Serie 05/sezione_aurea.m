function [xv, k, errv] = sezione_aurea(f, a, b, tol, nmax)

% [xv, k, errv] = sezione_aurea(f, a, b, tol, nmax)
% approssima il minimo di una funzione scalare f: [a,b] −> R tramite il
% metodo della sezione aurea
% IN
%   - f: funzione obiettivo
%   - a: estremo inferiore dell'intervallo
%   - b: estremo superiore dell'intervallo
%   - tol: tolleranza criterio d'arresto basato su stima dell'errore
%   - nmax: numero massimo di iterazioni ammesse
% OUT
%   - xv: vettore delle approssimazioni x_k del minimo
%   - k: numero di iterazioni eseguite
%   - errv: vettore delle stime dell'errore e

k = 0;

% salvo il rapporto aureo
phi = (1 + sqrt(5))/2;
% calcolo la prima stima
x = (a + b)/2;
% inizializzo il vettore delle soluzioni
xv = x;

% calcolo la stima iniziale dell'errore
err = (b - a)/2;
% inizializzo il vettore degli errori
errv = err;

% finché siamo sotto al numero di iterazioni massimo e se la stima
% dell'errore è maggiore della tolleranza
while k < nmax && err > tol
    k = k + 1;
    % calcolo c e d
    c = a + (b - a)/(phi + 1);
    d = a + (b - a)/phi;

    % verifico se siamo nell'intervallo di SX o di DX
    if (f(c) > f(d))
        a = c;
    else
        b = d;
    end

    % aggiorno la stima e la salvo nel vettore
    x = (a + b)/2;
    xv = [xv x];

    % aggiorno la stima dell'errore e salvo nel vettore
    err = (b - a)/2;
    errv = [errv err];
end

% printo le informazioni di convergenza o meno
if err <= tol
    fprintf("Il metodo della sezione aurea converge in %d iterazioni \n al punto di minimo x = %f \n", k, x);
else
    fprintf("Il metodo della sezione aurea non converge in %d iterazioni. \n", k)
end

end
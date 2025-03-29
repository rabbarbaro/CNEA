function [xvect, it] = bisez(a, b, toll, fun)

% [xvect, it] = bisez(a, b, toll, fun)
% metodo di bisezione per la risoluzione approssimata dell'equazione non
% lineare f(x) = 0 (ricerca dello zero)
% IN
%   - a, b: estremi dell'intervallo di ricerca
%   - toll: tolleranza sul test d'arresto (residuo)
%   - fun: funzione inline/anonymus
% OUT
% - xvect: vettore (colonna) con tutte le iterate calcolate (l'ultima è la
%       soluzione)
% - it: iterazioni effettuate

%% verifica applicabilità

if fun(a) * fun(b) > 0
    error ("La funzione deve avere segno diverso nei due estremi");
end

%% inizializzazione

xvect = [];
err = toll + 1;
nmax = ceil(log2((b-a)/toll) - 1);
it = -1;

%% algoritmo

while it < nmax && err > toll
    it = it + 1;

    x = (a + b) / 2;
    fun_c = fun(x);

    if fun_c == 0
        err = 0;
    else
        err = abs(fun_c);
    end

    xvect = [xvect; x];

    if fun(a) * fun_c < 0
        b = x;
    else
        a = x;
    end
end

if (it==nmax)
    fprintf("Massimo numero di iterazioni raggiunto. Errore sul residuo: %-6.4e \n", err);
else
    fprintf("x_%d soddisfa la tolleranza sul residuo \n", it);
end

fprintf("Zero calcolato: %-12.8f \n", xvect(end));

end
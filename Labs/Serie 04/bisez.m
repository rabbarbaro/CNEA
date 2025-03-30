function [xvect, it] = bisez(a, b, toll, fun)

% [xvect, it] = bisez(a, b, toll, fun)
% Metodo di bisezione per la risoluzione approssimata dell'equazione non
% lineare f(x) = 0 (ricerca dello zero)
% IN
%   - a, b: estremi dell'intervallo di ricerca
%   - toll: tolleranza sul test d'arresto (residuo)
%   - fun: funzione inline/anonymus
% OUT
%   - xvect: vettore (colonna) con tutte le iterate calcolate (l'ultima è
%       la soluzione)
% - it: iterazioni effettuate

%% verifica applicabilità

if fun(a) * fun(b) > 0
    error ("La funzione deve avere segno diverso nei due estremi");
end

%% inizializzazione

% inizializziamo il vettore delle iterate
xvect = [];
% errore arbitrario maggiore della tolleranza
err = toll + 1;

% calcoliamo dal risultato teorico il numero massimo di iterazioni in cui
% sicuramente convergeremo
nmax = ceil(log2((b-a)/toll) - 1);

% inizializziamo le iterazioni a -1, così la prima iterazione sarà in
% realtà il calcolo di x_0 (per non ripetere codice)
it = -1;

%% algoritmo

% finché siamo sotto al numero di iterazioni massimo e l'errore sul residuo
% è maggiore della tolleranza
while it < nmax && err > toll
    it = it + 1;
    % calcolo la nuova stima dello zero e valuto la funzione in quel punto
    x = (a + b) / 2;
    fun_c = fun(x);
    % se ho trovato lo zero pongo l'errore a 0 (uscita dal ciclo),
    % altrimenti calcolo l'errore sul residuo
    if fun_c == 0
        err = 0;
    else
        err = abs(fun_c);
    end
    
    % "appendo" la stima trovata al vettore
    xvect = [xvect; x];

    % se lo zero è all'interno del primo intervallo (a, x) lo seleziono,
    % altrimenti seleziono l'altro
    if fun(a) * fun_c < 0
        b = x;
    else
        a = x;
    end
end

% diamo informazioni nel terminale sul motivo di arresto e sul residuo
if (it==nmax)
    fprintf("Massimo numero di iterazioni raggiunto. Errore sul residuo: %-6.4e \n", err);
else
    fprintf("x_%d soddisfa la tolleranza sul residuo \n", it);
end

% printiamo l'ultima iterazione
fprintf("Zero calcolato: %-12.8f \n", xvect(end));

end
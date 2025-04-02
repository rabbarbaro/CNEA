function [xvect, it] = newton(x0, nmax, toll, fun, dfun, mol)

% [xvect, it] = newton(x0, nmax, toll, fun, dfun, (opt) mol)
% Metodo di Newton per la ricerca degli zeri della funzione fun, con
% criterio d'arresto basato sul controllo della differenza tra due iterate
% successive
% IN
%   - x0: stima iniziale
%   - nmax: numero massimo di iterazioni
%   - toll: tolleranza sul criterio d'arresto
%   - fun: funzione inline/anonymus
%   - dfun: derivata di fun come funzione inline/anonymus
%   - (opt) mol: opzionale, molteplicità dello zero per applicare il metodo
%       di Newton modificato
%       - se non assegnato viene imposto a 1 (metodo di Newton normale)
% OUT
%   - xvect: vettore (colonna) con tutte le iterate calcolate (l'ultima è
%       la soluzione)
%   - it: iterazioni effettuate

%% scelta newton o newton modificato

% se non viene passata la molteplicità la mettiamo a 1
if nargin == 5
    mol = 1;
end

%% inizializzazione

% il primo elemento nella soluzione è la stima iniziale x_0
xvect = x0;
% inzializziamo le iterate a 0 e un errore arbitrario maggiore della
% tolleranza per entrare nel ciclo
it = 0;
err = toll + 1;

%% algoritmo

% finché siamo l'errore è maggiore della tolleranza e non abbiamo raggiunto
% il numero massimo di iterazioni
while err > toll && it < nmax
    % se la derivata si annulla arrestiamo, dovremmo fare una divisione per
    % zero
    if dfun(x0) == 0
        error(' Arresto per azzeramento di dfun');
    else
        % calcolo la nuova iterata seguendo l'algoritmo
        x = x0 - mol * fun(x0) / dfun(x0);
        % calcolo la differenza fra iterate successive
        err = abs(x - x0);
        % aggiorno il vettore delle iterate
        xvect = [xvect; x];
        % salvo l'iterata corrente per l'iterazione successiva
        x0 = x;
        it = it + 1;
    end
end

if (it < nmax)
    fprintf("Convergenza al passo %d \n", it);
else
    fprintf("Raggiunto il numero massimo di passi %d \n", it);
end
fprintf("Radice calcolata: %-12.8f \n", xvect(end));

end
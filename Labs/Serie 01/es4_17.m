clear
clc
close all

n = 10;
v = randi(1e3, 1, n);

% chiamo la funzione + verifica corretto risultato
sorted = badSort(v);
% isequal(sorted, sort(v, 'descend'))

function v = badSort(v)
    n = length(v);
    % per ogni posizione i:
    % salvo il primo valore in una variabile temporanea
    % vedo qual è il massimo dei valori contentuti nelle posizioni da i:n
    % metto quel valore in posizione i e ci sostituisco la temporanea
    % aumento la i
    for ii = 1:n
        temp = v(ii);
        % se chiedo a max() due argument mi dà valore e posizione
        [maxVal, maxPos] = max(v(ii:n));
        v(ii) = maxVal;
        v(maxPos + ii - 1) = temp;
    end
end
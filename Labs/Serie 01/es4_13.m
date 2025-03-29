clear
clc
close all

% definisco la penultima cifra della matricola e la dimensione della
% tabella
penultimaMatricola = 8;
size = 10;

tabella = matricola2tabella(penultimaMatricola, size);

% funzione che cicla tutte le posizioni ed esegue i calcoli
function tab = matricola2tabella(k, n)
    tab = zeros(n);
    for ii = 1:n
        for jj = 1:n
            tab(ii, jj) = 2*ii*jj + (k + 1);
        end
    end
end
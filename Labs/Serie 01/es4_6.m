clear
clc
close all

% definisco la dimensione del vettore e prealloco a per efficienza del
% codice (sa la dimensione)
n = 10;
a = zeros(1, n);

% seguo i suggerimenti del testo per creare il vettore
for k = 1:n
    if k == 2 || k == 6
        a(k) = 1/k;
    else
        a(k) = 1 / ((k - 2) * (k - 6));
    end
end
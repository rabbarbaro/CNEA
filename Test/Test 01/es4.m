clear
clc
close all

% definisco la dimensione e prealloco la matrice
n = 100;
A = zeros(n);
% definisco c_k
k_max = n-1;
c_k = n*ones(1, 2*n - 1) - abs([-k_max:k_max]);

% si può fare anche così
c = n-(0:k_max);
toeplitz(c);

for ii = 1:n
    for jj = 1:n
        % devo aggiungere 100 all'indice, dato che c_-99 diventa c_1
        A(ii, jj) = c_k(ii-jj + 100);
    end
end

s1 = sum(A(n, :));
s2 = sum(diag(A(1:n, n:-1:1)));
% trattando la matrice come un vettore con le colonne in fila
% A(n:n-1:end-1) dà direttamente l'antidiagonale
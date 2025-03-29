clear
clc
close all

% costruisco la matrice di Hilbert 5x5 a mano
% nei cicli for non uso i, che è l'unità immaginaria
n = 5;
% prealloco A per efficienza del codice (sa la dimensione)
A = zeros(n);
for a = 1:n
    for b = 1:n
        A(a,b) = 1/(a + b - 1);
    end
end

% verifico la correttezza con il comando nativo di Matlab
% isequal(hilb(n), A);
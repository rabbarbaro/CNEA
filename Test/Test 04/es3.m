clear
clc
close all

n = 4;
A = 3*diag(ones(n, 1)) - diag(ones(n-1, 1), -1) + 2*diag(ones(n-1, 1), 1);
b = ones(n, 1);
x0 = b;
[x, k] = gs(A, b, x0, 1e-12, 5)

%%

A = [3 -1
    -1 2];
beta = 5/4;
P = [beta, -1
    0 2];

eig(eye(2) - P\A)

%%

A = hilb(4);
d = (cond(A)-1) / (cond(A)+1);

k_min = log(1e-3)/log(d)

%%

n = 100;
A = 8.1*diag(ones(n, 1)) - 3*diag(ones(n-1, 1), -1) - 3*diag(ones(n-1, 1), 1) ...
    - diag(ones(n-2, 1), -2) - diag(ones(n-2, 1), 2);

betav = 2:.001:2.5;

% inizializzo il vettore con i d per ogni beta
dv = [];

% scorrendo tutti i beta
for beta = betav
    % genero la matrice precondizionatore P
    P = diag(beta * ones(n, 1)) - diag(ones(n-1, 1), 1) ...
    - diag(ones(n-1, 1), -1);
    % calcolo il numero di condizionamento di P^-1 * A
    condinvPA = cond(P\A);
    % calcolo il valore d e lo "appendo" al vettore
    d = (condinvPA - 1 ) / (condinvPA + 1);
    dv = [dv, d];
end

% plotto per visualizzare meglio i risultati
plot(betav, dv)
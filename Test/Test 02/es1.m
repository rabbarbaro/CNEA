clear
clc
close all

% scelgo il valore e inizializzo la somma a 1
% metto in una matrice 2xiterazioni il valore della somma e l'iterazione
x = 2;
S(1, 1) = 1;
toll = 1e-1;
% inizializzo l'incremento a valore arbitrario > di toll
incr = 1;

% l'iterazione effettiva sarà k-1
k = 1;
% finché l'incremento non è minore della tolleranza
while incr > toll
    % salvo il numero dell'iterazione
    S(2, k+1) = k;
    % applico la somma come da testo
    S(1, k+1) = S(1, k) + x^k / factorial(k);
    % calcolo l'incremento
    incr = S(1, k+1) - S(1, k);
    k = k+1;
end

% in questo caso conosco la soluzione esatta, la uso per verificare
% S_2(1, 1) = 1;
% k = 1;
% while exp(x)-S_2(1, k) > toll
%     S_2(2, k+1) = k;
%     S_2(1, k+1) = S_2(1, k) + x^k / factorial(k);
%     incr = S_2(1, k+1) - S_2(1, k);
%     k = k+1;
% end
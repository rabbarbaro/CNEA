clear
clc
close all

K = 8;

vettore = kToVec(K);

% la funzione contiene l'esercizio 4_7
function vec = kToVec(K)
    vec = zeros(1, K+1);
    for k = 0:K
        vec(k + 1) = (2 * k + 1)^2;
    end
end
clear
clc
close all

% seleziono n_0 arbitrario
n_0 = 7;

% non è necessario dichiarare questa variabile n_k, una volta inizializzata
% la successione v(1) = n_0 posso direttamente iterare da lì
n_k = n_0;
% devo inizializzare il primo valore della successione
v(1) = n_0;

% applico l'algoritmo dell'esercizio
while n_k ~= 1
    if mod(n_k, 2) ~= 0
        n_k = 3 * n_k + 1;
    else
        n_k = n_k/2;
    end
    % uso uno dei possibili modi per aggiungere un elemento a un vettore
    % alternativamente esiste v = [v n_k]
    v(end+1) = n_k;
end
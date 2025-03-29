clear
clc
close all

B = [10 -1  1  0
      1  1 -1  3
      2  0  2 -1
      3  0  1  5];

%% 1

% costruiamo i cerchi di Gershgorin per avere una stima visiva della
% posizione degli autovalori
gershcircles(B)
% subito possiamo vedere che il massimo valore possibile del modulo
% degli autovalori è 12

%% 2

eig(B)
% gli autovalori di modulo minimo sono complessi coniugati, hanno lo stesso
% modulo (non distinto)

%% 3

% definisco (educated guess, non data da testo) tolleranza, numero massimo
% di iterazioni e vettore di partenza
tol = 1e-6;
nmax = 1000;
x0 = ones(4, 1);

% chiamo la funzione invpower
[lambda, x, iter] = invpower(B, tol, nmax, x0);
% in effetti il metodo non converge! (i due autovalori di modulo minimo
% hanno lo stesso modulo, sono complessi coniugati)

%% 4

% creiamo la funzione che implemennta il metodo delle potenze inverse con
% shift
% vedi file invpowershift.m

% per trovare i due autovalori con modulo minimo mi sposto dall'asse reale
% di +- 1i
s1 = 1i;
[lambda1, x1, iter1] = invpowershift(B, s1, tol, nmax, x0);
s2 = -1i;
[lambda2, x2, iter2] = invpowershift(B, s2, tol, nmax, x0);

% stampo i valori trovati
fprintf("autovalori c.c. di modulo minimo: \n%f + %fi \n%f + %fi \n", ...
    real(lambda1), imag(lambda1), real(lambda2), imag(lambda2))

% dato che dai cerchi di gershgorin troviamo che il massimo valore
% possibile del modulo degli autovalori è 12
shiftmax = 12;
[lambdamax1, xmax1, itermax1] = invpowershift(B, shiftmax, tol, nmax, x0);

fprintf("con invpowershift troviamo l'autovalore di modulo massimo in %d iterazioni: \n%f + %fi \n", ...
    itermax1, real(lambdamax1), imag(lambdamax1))

% confrontiamo con il metodo delle potenze
[lambdamax2, xmax2, itermax2] = eigpower(B, tol, nmax, x0);

fprintf("con eigpower troviamo l'autovalore di modulo massimo in %d iterazioni: \n%f + %fi \n", ...
    itermax2, real(lambdamax2), imag(lambdamax2))
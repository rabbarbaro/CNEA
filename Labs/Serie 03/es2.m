clear
clc
close all

%% 1

% creiamo la funzione che implemennta il metodo delle potenze inverse
% vedi file invpower.m

%% 2

% costruiamo la matrice A come da testo
A = toeplitz(1: 4);

% definiamo tolleranza, iterazioni massime e stima iniziale come da testo
tol = 1e-6;
nmax = 100;
x0 = [1:4]';

% chiamiamo il metodo delle potenze inverse
[lambda, x, iter] = invpower(A, tol, nmax, x0);

% stampiamo i risultati
if iter ~= nmax
     fprintf("invpower converge all'autovalore %f in %d iterazioni \n", ...
         lambda, iter);
else
     fprintf('invpower non converge in %d iterazioni. \n', iter)
end

% modifichiamo ora la stima iniziale
x0_1 = ones(4, 1);

% chiamiamo il metodo delle potenze inverse
[lambda_1, x_1, iter_1] = invpower(A, tol, nmax, x0_1);

% stampiamo i risultati
if iter ~= nmax
     fprintf("invpower converge all'autovalore %f in %d iterazioni \n", ...
         lambda_1, iter_1);
else
     fprintf('invpower non converge in %d iterazioni. \n', iter_1)
end

% converge ad un autovalore, ma non è quello di modulo minimo!
% la stima iniziale è perpendicolare all'autovalore che ci interessa
% verifica:
% [a, b] = eig(A);
% a(:, 3)' * x0_1;

%% 3

% per forzare la convergenza all'autovalore minimo nel caso in cui la stima
% iniziale fosse perpendicolare all'autovettore di interesse continuiamo a
% iterare senza fermarci una volta raggiunta una tolleranza (le variazioni
% indotte dalle approssimazioni numeriche causano uno spostamento
% dell'autovettore in una direzione non perpendicolare all'autovettore, da
% cui poi possiamo convergere all'autovalore minimo)

% vedi funzione invpower_mod.m
[lambda_2, x_2, min_eigenvec_comp] = invpower_mod(A, nmax, x0_1);
% quando min_eigenvec_comp = 1 abbiamo raggiunto convergenza all'autovalore
% minimo (corretto)

% stampiamo i risultati
if iter ~= nmax
     fprintf("invpower_mod converge all'autovalore %f\n", lambda_2);
else
     fprintf('invpower non converge in %d iterazioni. \n', iter_1)
end

subplot(1,2,1)
plot(min_eigenvec_comp)
subplot(1,2,2)
semilogy(min_eigenvec_comp)
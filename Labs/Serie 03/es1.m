clear
clc
close all

%% 1

% si definisce ranking (importanza relativa) r(p) della pagina p la somma
% di tutti i ranking r(q) di tutte le pagine q che puntano p
% r(p) = Sum_q->p (r(q)/#q))
% con #q numero di link presenti nella pagina q

% costruiamo la matrice di connessione A come dal testo
n_A = 100;

% 1. poniamo un uno o uno zero per ogni elemento (A)_ij che connette o meno
% l'elemento i-esimo al j-esimo (random)
A = randi([0 1], n_A);

% 2. l'elemento (A)_ij è o zero o 1/(numero di link presenti nella pagina),
% rappresentando la probabilità che si scelga a caso quel link (essendo un
% vettore che rappresenta una distribuzione di probabilità la somma sulle
% colonne deve fare sempre 1)
% (il comando sum somma sulle colonne)
s_A = sum(A);
for ii = 1:size(A, 1)
    A(ii, :) = A(ii, :) ./ s_A;
end

% il ranking delle pagine web p_i è rappresentato dal vettore colonna
% PageRank r = [r_1 r_2 ... r_n]'
% risolvere l'equazione iniziale equivale al problema r = Ar

% problema agli autovalori!
% il PageRank è quindi l'autovettore corrispondente all'autovalore 1 del
% problema agli autovalori associato.
% si puo dimostrare che se lambda_i sono gli autovalori di A allora
% abs(lambda_i) <= 1 e lambda_1 ha molteplicita uno (se ogni pagina
% ha almeno un link)

%% 2

% costruiamo la matrice come da testo
n_B = 5;
B = [0 0 0 1 0
    1 0 0 0 0
    0 1 0 0 0
    0 1 0 0 1
    1 1 1 1 0];

% normalizziamo le colonne
s_B = sum(B);
for ii = 1:size(B, 1)
    B(ii, :) = B(ii, :) ./ s_B;
end

%% 3

% creiamo la funzione che implemennta il metodo delle potenze
% vedi file eigpower.m

%% 4

% definiamo tolleranza, iterazioni massime e stima iniziale come da testo
tol = 1e-6;
nmax = 100;
x0_A = 1/n_A * ones(n_A, 1);

% chiamiamo il metodo delle potenze
[lambda_A, x_A, iter_A] = eigpower(A, tol, nmax, x0_A);

% stampiamo i risultati
if iter_A ~= nmax
     fprintf("eigpower converge all'autovalore %f in %d iterazioni \n", ...
         lambda_A, iter_A);
else
     fprintf('eigpower non converge in %d iterazioni. \n', iter_A)
end

% facciamo la stessa operazione per la matrice B

x0_B =  1/n_B * ones(n_B, 1);
[lambda_B, x_B, iter_B] = eigpower(B, tol, nmax, x0_B);

if iter_B ~= nmax
     fprintf("eigpower converge all'autovalore %f in %d iterazioni \n", ...
         lambda_B, iter_B);
else
     fprintf('eigpower non converge in %d iterazioni. \n', iter_B)
end
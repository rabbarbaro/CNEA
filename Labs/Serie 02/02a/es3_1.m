clear
clc
close all

% abbiamo N molle ed N+1 variabili:
% - x_0 = 0 (1)
% - x_i per i = 1:N-1 (N-1) posizione dell'anello i-esimo
% - x_N = L (1)
% ovvero N-1 incognite (dimensione del problema)
% scriviamo N-1 equazioni, ovvero i bilanci delle forze per ciascun anello
% (siamo in configurazione di equilibrio, quindi ci sarà la forza
% esercitata dalla molla a SX che bilancia la forza esercitata dalla molla
% a DX di ogni anello)
% - F_i + F_i+1 = 0 per i = 1:N-1 (N-1)
% la forza esercitata da una molla è proporzionale al suo allungamento 
% rispetto alla sua lunghezza in una configurazione a riposo (F = K * dL)
% Utilizzando le variabili x_i, considerando nulla la lunghezza a riposo,
% l’equazione di equilibrio perl'i-esimo nello è:
% -K(x_i - x_i-1) + K(x_i+1 - x_i) = 0
% K*x_i-1 - 2K*x_i + K*x_i+1 = 0

% dato che:
% - x_0 = 0 (1)
% - x_N = L (1)
% possiamo scrivere il sistema complessivo come:
% A = K * [-2  1  0  ...   0  0
%           1 -2  1  ...   0  0
%           0  1 -2  ...   0  0
%           0  0  0  ...  -2  1
%           0  0  0  ...   1 -2]
% b = [0 0 ... 0 -K*L]'
% dato che l'allungamento della prima molla sarà x_i-0 e quello dell'utlima
% molla sarà L-X_N-1

% inserisco i dati
K = 100;
L = 20;
N = 20;

% costruisco A
n = N-1;
A = diag(-2*(ones(n, 1))) + diag(ones(n-1, 1), -1) + diag(ones(n-1, 1), 1);
A = K*A;

% definisco il vettore f
f = zeros(n,1);
f(end) = -K * L;

% chiamo le matrici LL ed UU per evitare equivoci con la lunghezza L
[LL, UU, x] = thomas(A, f);

% verifichiamo che il risltato trovato con la nostra implementazione di
% Thomas sia uguale al risultato che trova Matlab
xx = A \ f;
norm(x - xx) / norm(x);

% se a ciascun anello viene esercitata una forza esterna F_i_ext
% l'equazione di equilibrio per ogni anello si modifica in:
% -K(x_i - x_i-1) + K(x_i+1 - x_i) + F_i_ext = 0
% K*x_i-1 - 2K*x_i + K*x_i+1 = -F_i_ext
% dunque la matrice A non subisce variazioni, mentre il vettore b:
% b = [-F_1_ext -F_2_ext ... -F_N-2_ext -F_N-1_ext-K*L]'
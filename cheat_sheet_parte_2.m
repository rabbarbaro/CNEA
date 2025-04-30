%% interpolazione polinomiale

% matlab salva i polinomi come array dei coefficienti
% p = [a_n, a_n-1, ..., a1, a0]

% per ottenere il polinomio interpolatore di lagrange passante per le n+1
% coppie di dati {x_i, y_i} con i = 0, ..., n
% x array dei nodi
% y array dei valori nei nodi
% n grado del polinomio (max numero dei punti - 1)
p = polyfit(x, y, n);

% per valutare il polnomio in uno o più punti
% z array dei punti in cui valutarlo
pz = polyval(p, z);

% n+1 nodi equispaziati
x_nod = linspace(a, b, n+1);
% oppure
h = (b - a)/n;
x_nod = a:h:b;

% n+1 nodi di Chebyshev
tv = -cos(pi * (0:n)/n);
x_nod_che = (b - a)/2 * tv + (a + b)/2;

% polinomio interpolatore di grado n su nodi equispaziati
x = linspace(a, b, 1000);       % 1000 numero arbitrario di punti
% creo nodi equispaziati
x_nod = linspace(a, b, n+1);
% valuto la funzione nei nodi
f_nod = f(x_nod);
% ottengo il polimonio di lagrange dai nodi ottenuti e lo valuto su [a, b]
P = polyfit(x_nod, f_nod, n);
P_dis = polyval(P, x);

% polinomio interpolatore di grado n su nodi di Chebyshev
x = linspace(a, b, 1000);       % 1000 numero arbitrario di punti
% creo un vettore con i parametri t_k per k = 0, ..., n
% di fatto sono gli n+1 nodi di Chebyshev nell'intervallo di
% riferimento [-1, 1]
tv = -cos(pi * (0:n)/n);
% mappo i nodi dall'intervallo di riferimento ad [a, b]
x_nod_che = (b - a)/2 * tv + (a + b)/2;
% valuto la funzione nei nodi
f_nod_che = f(x_nod_che);
% ottengo il polimonio di lagrange dai nodi ottenuti e lo valuto
% sull'intervallo definito prima
P_che = polyfit(x_nod_che, f_nod_che, n);
P_dis_che = polyval(P_che, x);

%% approssimazione polinomiale nel senso dei minimi quadrati

% dai n+1 coppie di dati troviamo l'approssimazione polinomiale nel senso
% dei minimi quadrati
% x array dei nodi
% y array dei valori nei nodi
% m grado del polinomio < n
p = polyfit(x, y, m);
% poi da valutare con polyval

%% interpolante lineare a tratti

% approssimazione polinomiale a tratti, non va valutata ma restituisce già
% il valore su tutti i punti di z
% x array dei nodi
% y array dei valori nei nodi
% z array dei punti su cui valutare
pz = interp1(x, y, z);

%% splines

% spline naturale cubica
% non va valutata ma restituisce già il valore su tutti i punti di z
% x array dei nodi
% y array dei valori nei nodi
% z array dei punti su cui valutare
pz = cubicspline(x, y, z);

% spline not-a-knot
%, non va valutata ma restituisce già il valore su tutti i punti di xq
% x array dei nodi
% y array dei valori nei nodi
% xq array dei punti su cui valutare
s = spline(x, y, xq);

%% da aggiornare
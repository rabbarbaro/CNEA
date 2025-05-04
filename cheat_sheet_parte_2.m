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

%% calcolo della base lagrangiana

...

% fino al grado 2:
%
% % vettore dei nodi
% x_v = [0, 0.5, 2];
% 
% % applico la definizione dei polinomi base di Lagrange
% Phi1 = @(x) (x - x_v(2)) ./ (x_v(1) - x_v(2)) .* (x - x_v(3)) ./ (x_v(1) - x_v(3));
% Phi2 = @(x) (x - x_v(1)) ./ (x_v(2) - x_v(1)) .* (x - x_v(3)) ./ (x_v(2) - x_v(3));
% Phi3 = @(x) (x - x_v(1)) ./ (x_v(3) - x_v(1)) .* (x - x_v(2)) ./ (x_v(3) - x_v(2));
% 
% % plotto la funzione e i polinomi base
% plot(x, f(x), 'LineWidth', 2)
% hold on
% plot(x, Phi1(x), x, Phi2(x), x, Phi3(x))
% plot(x_v, 0*x_v, '*')
% plot(x_v, 1, 'k*')
% 
% % costruisco il polinomio interpolante di lagrange e plotto
% p = @(x) f(x_v(1)) .* Phi1(x) + f(x_v(2)) .* Phi2(x) + f(x_v(3)) .* Phi3(x);

% per il grado 3 vedi Test 08\es3.m

%% calcolo e stima costante di lebesgue

...

%% errore di interpolazione

err = abs(f_ex - f_interp);
err_max = max(err) = max(abs(f_ex - f_interp));

% dalla teoria la stima a priori dell'errore di interpolazione con lagrange su n+1 nodi equispaziati è:
%   - H = (b-a) / n
%   - E_H <= H^2 / 8 * max|f''(x)|, x in [a,b]
d2f = ...; % derivata seconda di f
max_d2f = max(abs(d2f(linspace(a, b, 1000)))); % massimo della derivata seconda di f in [a,b]
H = (b - a) / 2; % passo di interpolazione
E_H = H^2 / 8 * max_d2f; % errore di interpolazione

%% minimi quadrati non lineari

% volendo minimizzare la somma dei quadrati delle differenze tra i dati e
% il modello, si ha:
%   sum_i=0^n (y_i - f_tilde(x_i; w_vec))^2

% min_(x_vec in R^m) Phi(x_vec)
% Phi(y_vec) = 1/2 ||F_vec(y_vec)||^2 = 1/2 * sum_i=0^n (f_i(x_i))^2
% con il gradiente e l'hessiano definiti come:
%   grad_Phi(y_vec) = J(y_vec)^T * F(y_vec)
%   Hess_Phi(y_vec) = J(y_vec)^T * J(y_vec) + Q(y_vec)
% usiamo questo gradiente e hessiano per calcolare il minimo usando
% l'hessiana semplificata eliminando Q abbiamo il metodo di Gauss-Newton

% implementiamo la funzione gauss_newton che calcola il minimo di una
% funzione non lineare usando il metodo di Gauss-Newton
% vedi gauss_newton.m

% definisco i dati come colonna
x = linspace(...)';
% calcoliamo il valore della funzione con i dati veri per avere il valore
% nei nodi, definito come colonna
y = f(x, w_true);

% definisco il vettore di funzioni f (coi vari dati in colonna)
f = @(x, w); % funzione di x (nodi) e w(1), w(2), ecc...
% definisco il vettore dei residui (coi vari dati in colonna)
r = @(x, y, w) y - f(x, w);
% definisco lo jacobiano (derivo i residui, è una matrice)
J = @(x, y, w) [, ...
                ];
% lo jacobiano devo definirlo per riga

% definisco la guess iniziale, la tolleranza e il massimo numero di
% iterazioni
w0 = [...]';
toll = 1e-6;
maxit = 1000;
% chiamiamo la funzione
[w_vect, it] = gauss_newton(x, y, r, J, w0, toll, maxit);
w_vect(:, end)
% plottiamo
plot(x, y, 'o', x, f(x, w_vect(:, end)));
% notiamo che abbiamo effettivamente ritrovato i coefficienti corretti,
% dato che la legge trovata si sovrappone ai dati

% COME ESEMPIO VEDI:
%   - Lab\Serie 06\es8.m
%   - Test\Test 08\es9.m

%% formule di interpolazione trigonometriche (DFT)

% vedi Test 08\es6.m per n pari

%% formule di integrazione semplici

% dalla teoria:
%   - PUNTO MEDIO:
%       - punti medi:                   xm_k = (x_(k-1) + x_k) / 2
%       - lunghezza sttointervallo:     h = (b - a) / M
%           - numero di intervalli:     M
%       - formula:                      I = h * sum_k=1^M ( f(xm_k) )
%   - TRAPEZI:
%       - nodi:                         x_k
%       - lunghezza sttointervallo:     h = (b - a) / M
%       - formula:                      I = (h/2) * (f(a) + f(b)) + h * sum_k=1^(M-1) ( f(x_k) )
%   - SIMPSON:
%       - nodi:                         x_k
%       - punti medi:                   xm_k = (x_(k-1) + x_k) / 2
%       - lunghezza sttointervallo:     h = (b - a) / M
%       - formula:                      I = (h/6) * sum_k=1^M ( f(x_(k-1)) + 4 * f(xm_k) + f(x_k) )

% vedi le funzioni pmedcomp.m, trapcomp.m e simpcomp.m

%% errori formule di integrazione composite

% dalla teoria sappiamo che gli errori possono essere stimati come:
% PM e TRAP:
%   err <= C * ((b - a)/N)^2 * |f''|_max <= toll
%   N = ceil( (b-a) * sqrt((|f''|_max * C / toll)))
% con:
%   - C = (b - a) / 24      per PM
%   - C = (b - a) / 12      per TRAP
% SIMP
%   err <= C * ((b - a)/N)^4 * |f''''|_max <= toll
%   N = ceil( (b-a) * ((|f''''|_max * C / toll))^(1/4))
% con:
%   - C = (b - a) / (16 * 180)

% per il calcolo del numero minimo di intervalli:
% % definisco la tolleranza
% toll = 1e-6;
% % calcolo la derivata seconda e quarta della funzione di integrazione
% d2fun = @(x) ...;
% d4fun = @(x) ...;
% % % definisco i valori di C per i metodi del punto medio, trapezi e Simpson
% C_pm = (b - a) / 24;
% C_trap = (b - a) / 12;
% C_simp = (b - a) / (16 * 180);
% % % calcolo il numero minimo di intervalli per ogni metodo di integrazione in
% % modo che l'errore sia minore della tolleranza
% N_pm = ceil( (b - a) * sqrt((C_pm * max(abs(d2fun(x_q)))) / toll));
% N_trap = ceil( (b - a) * sqrt((C_trap * max(abs(d2fun(x_q)))) / toll));
% N_simp = ceil( (b - a) * ((C_simp * max(abs(d4fun(x_q)))) / toll)^(1/4));

%% formula di integrazione di Gauss-Legendre

% FORMULA A DUE PUNTI:
% dalla teoria, sull'intervallo di riferimento [-1, 1] e per la formula a
% due punti:
%   - nodi:         y_0_ = -1/sqrt(3);      y_1_ = +1/sqrt(3)
%   - pesi:         w_0_ = 1;               w_1_ = 1
% passando all'intervallo generico [a, b], sempre per due punti:
%   - cambio di variabile:       phi(y) = (a+b)/2 + (b-a)/2 * y
%   - nodi:                     y_0 = (a+b)/2 + (b-a)/2 * y_0_;
%                               y_1 = (a+b)/2 + (b-a)/2 * y_1_
%   - pesi:                     w_0 = (b-a)/2 * w_0_;
%                               w_1 = (b-a)/2 * w_1_
% dunque:
%   - nodi:         y_0 = (a+b)/2 - (b-a)/2 * 1/sqrt(3);
%                   y_1 = (a+b)/2 + (b-a)/2 * 1/sqrt(3)
%   - pesi:         w_0 = (b-a)/2;          w_1 = (b-a)/2
% in generale vale che:
%   I = sum_k=0^N ( w_k * f(y_k) )
% per quella a due punti
%   I = sum_k=0^1 ( w_k * f(y_k) )

% vedi la funzione gausscomp.m

% PER TROVARE NODI E PESI
% useremo la funzione zplege, che trova nodi e pesi di Gauss-Legendre
% definisco l'ordine della formula e il numero di nodi
%   - n indica l’ordine della formula
%   - (n + 1 sono il numero di nodi di quadratura)
[x, w] = zplege(n, a, b);
I = sum(w .* p1(x));

%% integrazione di clenshaw-curtis

% per n pari vedi Test 09\es4.m
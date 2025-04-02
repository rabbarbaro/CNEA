clear
clc
close all

% definisco la funzione come da testo
f = @(x) cos(2*x).^2 - x.^2;

%% 1

% consideriamo questo intervallo arbitrario
x = linspace(-1, 1, 1000);

% plottiamo
plot(x, f(x))
grid on
hold on
plot(x, 0*x, '--')

% la funzione è pari, basta che cerco lo zero positvo e trovo anche lo zero
% negativo automaticamente

%% 2

% consideriamo la funzione di iterazione
% phi(x) = x + A*f(x)

% per verificare in modo teorico l'intervallo di valori di A tali per cui
% il metodo delle iterazioni di punto fisso converge allo zero di (x)
% dobbiamo verificare la condizione
% |phi'(alpha)| < 1

% calcolando
% phi'(alpha) = 1 + A*f'(alpha)
% |1 + A*f'(alpha)| < 1
% -2 < A*f'(alpha) < 0
% f'(allpha) è negativo per l'alpha positivo, devo cambiare verso della
% disequazione quando divido

% 0 < A < 2/f'(alpha)

%% 3

% definisco tolleranza, massimo numero di iterazioni e intervallo
tol = 1e-10;
nmax = 1000;
a = -1;
b = 1;

% guess iniziale
x0 = 0.1;

% calcolo l'iterazione di punto fisso A = 0.1
A1 = 0.1;
phi1 = @(x) x + A1*f(x);
[succ1, it1] = ptofis(x0, phi1, nmax, tol, a, b);
alpha = succ1(end);

% ora che ho trovato una buona apporssimazione di alpha posso calcolare il
% valore della derivata di f in alpha
df = @(x) -4*cos(2*x).*sin(2*x) - 2*x;
% ottengo dunque l'estremo superiore per la scelta di alpha
Asup = -2/df(alpha)

% A2 < Asup
A2 = 0.6;
phi2 = @(x) x + A2*f(x);
[succ2, it2] = ptofis(x0, phi2, nmax, tol, a, b);

% A3 > sup
A3 = 0.75;
phi3 = @(x) x + A3*f(x);
[succ3, it3] = ptofis(x0, phi3, nmax, tol, a, b);
% come previsto dalla teoria NON converge

% A4 = A2, cambio x0
A4 = 0.6;
x0 = 2;
phi4 = @(x) x + A4*f(x);
[succ4, it4] = ptofis(x0, phi4, nmax, tol, a, b);
% per questo tipo di funzioni phi la scelta di x0 è meno rilevante rispetto
% alla velocità di convergenza rispetto al parametro A

%% 4

% usando stimap stimiamo l'ordine di convergenza e il fattore di
% convergenza del metodo di punto fisso
[p1, c1] = stimap(succ1);
[p2, c2] = stimap(succ2);
[p3, c3] = stimap(succ3);
[p4, c4] = stimap(succ4);
% quelle che convergono convergono tutte con ordine 1, e notiamo che il
% fattore A ha un effetto molto maggiore sul fattore di riduzione rispetto
% alla scelta di x0
% quello che non converge presenta un NaN

%% 5

% per avere convergenza del secondo ordine è necessario che
% phi'(alpha) = 0

% phi'(alpha) = 1 + A*f'(alpha) = 0
% Aopt = -1/f'(alpha)
% (ovvero Asup/2)

% calcoliamo A ottimo trovato dalla teoria e poi lanciamo il metodo
Aopt = -1/df(alpha);
phi_opt = @(x) x + Aopt*f(x);
[succ_opt, it_opt] = ptofis(0.1, phi_opt, nmax, tol, a, b);

% ci asoettiamo unordine di convergenza p = 2
[p_opt, c_opt] = stimap(succ_opt);

%% 6

% implementiamo il metodo di newton come iterazione di punto fisso
phi_N = @(x) x - f(x)./df(x);
[succ_N, it_N] = ptofis(0.1, phi_N, nmax, tol, a, b);
[p_N, c_N] = stimap(succ_N);
% effettivamente essendo uno zero semplice troviamo che im metodo di newton
% converga con ordine 2

% creiamo un vettore di possibili guess iniziali da 0.01 a 1 (non zero!) e
% lanciamo il metodo di newton per ognuna
vett_x0 = 0.01:.01:1;
vett_sol = [];
for x0 = vett_x0
    [succ_n, it_n] = ptofis(x0, phi_N, nmax, tol, a, b);
    vett_sol = [vett_sol succ_n(end)];
end

% plottiamo dove converge il metodo VS x0
plot(vett_x0, vett_sol)
% troviamo sperimentalmente che per convergere allo zero positivo
% 0.12 <= x0 <= 0.88
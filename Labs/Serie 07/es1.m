clear
clc
close all

% definiamo la funzione che descrive la sezione del solido nel piano
f = @(x) cosh(x - 0.5);

% definiamo l'intervallo di integrazione
a = 0;
b = 0.8;
x_q = linspace(a, b, 100);

% definiamo il volume esatto del solido
% V = pi * int_a^b f(x)^2 dx
V_ex = pi * ( (sinh(1) + sinh(3/5)) / 4 + 2/5);

%% 1

% plottiamo la funzione f(x)
plot(x_q, f(x_q));

%% 2, 3, 4

% scriviamo le funzioni di integrazione per i metodi del punto medio,
% trapezi e Simpson

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

%% 5

% inizializzo i vettori per i risultati delle integrazioni
I_pm = [];
I_t = [];
I_s = [];

% definisco il numero di intervalli
N = 1:20;

% definisco la funzione di cui calcolare l'integrale di volume
V_N = @(x) pi * f(x).^2;

% calcolo il volume del solido per ogni metodo di integrazione e per ogni
% numero di intervalli
for n = N
    I_pm = [I_pm, pmedcomp(a, b, n, V_N)];
    I_t = [I_t, trapcomp(a, b, n, V_N)];
    I_s = [I_s, simpcomp(a, b, n, V_N)];
end

% plottiamo i risultati delle integrazioni
figure
plot(N, I_pm, 'r', N, I_t, 'g', N, I_s, 'b', N, V_ex * ones(1, length(N)), 'k--')
legend('Punto medio', 'Trapezi', 'Simpson', 'Valore esatto')

%% 6

% calcoliamo gli errori relativi tra i risultati delle integrazioni e il
% valore esatto del volume del solido
err = [abs(I_pm - V_ex); abs(I_t - V_ex); abs(I_s - V_ex)];

% calcoliamo la lunghezza degli intervalli
H = (b - a)./N;

% plottiamo gli errori relativi in scala logaritmica in funzione della
% lunghezza degli intervalli
figure
loglog(H, err)
hold on
loglog(H, H.^2, '--', H, H.^4, '--')
legend('Err. Punto medio', 'Err. Trapezi', 'Err. Simpson', 'H^2', 'H^4')

%% 7

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

% definisco la tolleranza
toll = 1e-5;

% calcolo la derivata seconda e quarta della funzione di integrazione
% V_N = pi * f(x)^2
d2fun = @(x) 2*pi*cosh(x-0.5).^2 + 2*pi*sinh(x-0.5).^2;
d4fun = @(x) 8*pi*cosh(x-0.5).^2 + 8*pi*sinh(x-0.5).^2;

% definisco i valori di C per i metodi del punto medio, trapezi e Simpson
C_pm = (b - a) / 24;
C_trap = (b - a) / 12;
C_simp = (b - a) / (16 * 180);

% calcolo il numero minimo di intervalli per ogni metodo di integrazione in
% modo che l'errore sia minore della tolleranza
N_pm = ceil( (b - a) * sqrt((C_pm * max(abs(d2fun(x_q)))) / toll));
N_trap = ceil( (b - a) * sqrt((C_trap * max(abs(d2fun(x_q)))) / toll));
N_simp = ceil( (b - a) * ((C_simp * max(abs(d4fun(x_q)))) / toll)^(1/4));

fprintf('Numero minimo di intervalli per il metodo del punto medio: %d\n', N_pm)
fprintf('Numero minimo di intervalli per il metodo dei trapezi: %d\n', N_trap)
fprintf('Numero minimo di intervalli per il metodo di Simpson: %d\n', N_simp)
fprintf('Tolleranza: %e\n', toll)

% calcolo il volume del solido con i metodi del punto medio, trapezi e
% Simpson con il numero minimo di intervalli trovato
V_pm_tol = pmedcomp(a, b, N_pm, V_N);
err_pm = abs(V_pm_tol - V_ex);
V_trap_tol = trapcomp(a, b, N_trap, V_N);
err_trap = abs(V_trap_tol - V_ex);
V_simp_tol = simpcomp(a, b, N_simp, V_N);
err_simp = abs(V_simp_tol - V_ex);

fprintf('Errore relativo con il metodo del punto medio: %e\n', err_pm)
fprintf('Errore relativo con il metodo dei trapezi: %e\n', err_trap)
fprintf('Errore relativo con il metodo di Simpson: %e\n', err_simp)
% gli errori sono tutti sotto alla tolleranza!

%% 8

% scriviamo la funzione di integrazione per il metodo di Gauss-Legendre

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

%% 9

% inizializzo il vettore per i risultati delle integrazioni
V_gl = [];

% per lo stesso numero di intervalli N usato per gli altri metodi
for n = N
    V_gl = [V_gl, gausscomp(a, b, n, V_N)];
end

% plottiamo i risultati delle integrazioni
figure
plot(N, V_gl, 'm', N, V_ex * ones(1, n), 'k--')
legend('Gauss-Legendre', 'Valore esatto')

% plottiamo gli errori relativi in scala logaritmica in funzione della
% lunghezza degli intervalli
figure
loglog(H, abs(V_gl - V_ex))
hold on
loglog(H, H.^2, '--', H, H.^4, '--')
legend('Gauss-Legendre', 'H^2', 'H^4')

% l'ordine di convergenza di gauss-legendre è 4 (l'errore decresce come
% H^4)
% confrontato con trapezio e simpson:
%   - trapezio: ordine di convergenza 2
%       - formula a due punti come GL ma con ordine di convergenza pià
%         basso
%   - simpson: ordine di convergenza 4
%       - ha lo stesso ordine di convergenza di GL, ma ha un punto in più
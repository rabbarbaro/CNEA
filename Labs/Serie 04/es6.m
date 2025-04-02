clear
clc
close all

% definisco la funzione e il suo zero
f = @(x) cos(pi*x) .* (x - 0.5);
df = @(x) -pi * sin(pi*x) .* (x - 0.5) + cos(pi*x);

alpha = 0.5;

%% 1

% data tolleranza, numero massimo di iterazioni e guess iniziale
toll = 1e-6;
nmax = 1000;
x0 = 0.9;

% chiamo newton
[xvect_N, it_N] = newton(x0, nmax, toll, f, df, 1);

% stimo p e c
[p_N, c_N] = stimap(xvect_N);

% dalla teoria
% lim ||x^(k+1)-alpha||/||x^(k)-alpha||^p = 1/p! * d(p)phi(alpha)

% dato che p = 1
% lim ||x^(k+1)-alpha||/||x^(k)-alpha|| = dphi(alpha)

% dalla teoria
% dphi_newton(alpha) = 1 - 1/m

% m = 1 / (1 - dphi_newton(alpha))
m = 1/(1 - c_N(end));

% applico newton modificato
[xvect_Nm, it_Nm] = newton(x0, nmax, toll, f, df, m);

%% 2

% il criterio d'arresto applicato a newton semplice è soddisfacente se m =
% 1, per newton modificato è sempre soddisfacente

%% 3

% implementiamo la funzione che esegue il metodo delle secanti
% vedi file secanti.m

% hardcodo la prima iterazione (normalmente calcolerei con le corde)
x1 = 0.7;
xStart = [x0, x1];
nmaxSec = 10;

% calcolo con le secanti
[xvect_S, it_S] = secanti(xStart, nmaxSec, toll, f);

[p_S, c_S] = stimap(xvect_S);
% ci aspettiamo una convergenza lineare perché lo zero è multiplo (sarebbe
% stato (1 + sqrt(5)) / 2

%% 4

mu = 2;

% definisco phi e la sua derivata
phi = @(x) x + mu/(2*pi) * cos(pi * x);
dphi = @(x) 1 - mu/2 * sin(pi * x);

% per verificare in modo teorico l'intervallo di valori di A tali per cui
% il metodo delle iterazioni di punto fisso converge allo zero di (x)
% dobbiamo verificare la condizione
% |phi'(alpha)| < 1

% calcolando
% phi'(alpha) = 1 - mu/2 * sin(pi * alpha)
% | 1 - mu/2 * sin(pi * alpha) | < 1
% -2 < - mu/2 * sin(pi * alpha) < 0
% 0 < mu/2 * sin(pi * alpha) < 2
% 0 < mu < 4

% verifica
[succ, it] = ptofis(x0, phi, nmax, toll, -1, 1);

% per avere convergenza almento ordine 2
% phi'(alpha) = 0
% 1 - mu/2 * sin(pi * alpha) = 0
% mu/2 * sin(pi * alpha) = 1
% mu = 2

[p_PF, c_PF] = stimap(succ);

%%

% 0 < phi'(alpha) < 1
% 0 < 1 - mu/2 * sin(pi * alpha) < 1
% -1 < - mu/2 * sin(pi * alpha) < 0
% 0 < mu/2 * sin(pi * alpha) < 1
% 0 < mu/2 < 1
% 0 < mu < 2

%%

mu = 1;

% ridefinisco phi e derivata con il nuovo valore di mu
phi = @(x) x + mu/(2*pi) * cos(pi * x);
dphi = @(x) 1 - mu/2 * sin(pi * x);

% 0 < mu/2 * sin(pi * x) < 2
% mu = 1
% 0 < 1/2 * sin(pi * x) < 2
% 0 < sin(pi * x) < 4
% sin(pi * x) > 0           0 < pi * x < pi             0 < x < 1
% sin(pi * x) < 4           sempre verificata

% per: -1/2 <= a < alpha < b <= 3/2
% phi(x) in [a,b] per ogni x in [a,b] dato che:
% phi(a) >= a, phi(b) <= b e phi'(x) > 0, per ogni x in [−1/2,3/2]

% sistema
% -1/2 < x < 3/2
% 0 < x < 1

% 0 < x < 1
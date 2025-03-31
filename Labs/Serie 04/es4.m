clear
clc
close all

f = @(x) cos(2*x).^2 - x.^2;
x = linspace(-1, 1, 1000);

%% 1

plot(x, f(x))
grid on
hold on
plot(x, 0*x, '--')

% funzione pari, basta che cerco quello positvo e trovo anche lo zero
% negativo

%% 2

% phi(x) = x + A*f(x)
% phi'(alpha) = 1 + A*f'(alpha)
% |phi'(alpha)| < 1
% |1 + A*f'(alpha)| < 1
% -2 < A*f'(alpha) < 0
% f'(allpha) è negativo per l'alpha positivo
%2 0 < A < 2/f'(alpha)

%% 3

x0 = 0.1;
tol = 1e-10;
nmax = 1000;
a = -1;
b = 1;

A1 = 0.1;
phi1 = @(x) x + A1*f(x);
[succ1, it1] = ptofis(x0, phi1, nmax, tol, a, b);
alpha = succ1(end);

df = @(x) -4*cos(2*x).*sin(2*x) - 2*x;
% estremo superiore per la scelta di alpha
Asup = -2/df(alpha);

A2 = 0.6;
phi2 = @(x) x + A2*f(x);
[succ2, it2] = ptofis(x0, phi2, nmax, tol, a, b);

A3 = 0.75;
phi3 = @(x) x + A3*f(x);
[succ3, it3] = ptofis(x0, phi3, nmax, tol, a, b);

A4 = 0.6;
% cambio x0
phi4 = @(x) x + A4*f(x);
[succ4, it4] = ptofis(2, phi4, nmax, tol, a, b);

%% 4

[p1, c1] = stimap(succ1);
[p2, c2] = stimap(succ2);
[p3, c3] = stimap(succ3);
[p4, c4] = stimap(succ4);

%% 5

% per avere convergenza del secondo ordine phi'(alpha) = 0
% phi'(alpha) = 1 + A*f'(alpha) = 0
% Aopt = -1/f'(alpha)

Aopt = -1/df(alpha);

phi_opt = @(x) x + Aopt*f(x);
[succ_opt, it_opt] = ptofis(0.1, phi_opt, nmax, tol, a, b);
[p_opt, c_opt] = stimap(succ_opt);

%% 6

phi_n = @(x) x - f(x)./df(x);
[succ_n, it_n] = ptofis(0.1, phi_n, nmax, tol, a, b);
[p_n, c_n] = stimap(succ_n);

vett_x0 = 0.01:.01:1;
vett_sol = [];
for x0 = vett_x0
    [succ_n, it_n] = ptofis(x0, phi_n, nmax, tol, a, b);
    vett_sol = [vett_sol succ_n(end)];
end

plot(vett_x0, vett_sol)

% x0 in [0.12, 0.88]
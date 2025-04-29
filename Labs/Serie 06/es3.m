clear
clc
close all

% definisco funzione e intervallo

f = @(x) -x.^3 + 3*x.^2 - 2;
a = 0;
b = 2;
x = linspace(a, b, 1000);

%% 1

% vettore dei nodi
x_v = [0, 0.5, 2];

% applico la definizione dei polinomi base di Lagrange
Phi1 = @(x) (x - x_v(2)) ./ (x_v(1) - x_v(2)) .* (x - x_v(3)) ./ (x_v(1) - x_v(3));
Phi2 = @(x) (x - x_v(1)) ./ (x_v(2) - x_v(1)) .* (x - x_v(3)) ./ (x_v(2) - x_v(3));
Phi3 = @(x) (x - x_v(1)) ./ (x_v(3) - x_v(1)) .* (x - x_v(2)) ./ (x_v(3) - x_v(2));

% plotto la funzione e i polinomi base
plot(x, f(x), 'LineWidth', 2)
hold on
plot(x, Phi1(x), x, Phi2(x), x, Phi3(x))
plot(x_v, 0*x_v, '*')
plot(x_v, 1, 'k*')

% costruisco il polinomio interpolante di lagrange e plotto
p = @(x) f(x_v(1)) .* Phi1(x) + f(x_v(2)) .* Phi2(x) + f(x_v(3)) .* Phi3(x);

plot(x, p(x), 'LineWidth', 2)

%% verifica

% uso la funzione built-in

P = polyfit(x_v, f(x_v), 2);
P_dis = polyval(P, x);

figure
plot(x, P_dis, 'LineWidth', 2)
hold on
plot(x, p(x))

%% 2

% vettore dei nodi equispaziati
x_v_eq = [0, 1, 2];

% applico la definizione dei polinomi base di Lagrange
Phi1_eq = @(x) (x - x_v_eq(2)) ./ (x_v_eq(1) - x_v_eq(2)) .* (x - x_v_eq(3)) ./ (x_v_eq(1) - x_v_eq(3));
Phi2_eq = @(x) (x - x_v_eq(1)) ./ (x_v_eq(2) - x_v_eq(1)) .* (x - x_v_eq(3)) ./ (x_v_eq(2) - x_v_eq(3));
Phi3_eq = @(x) (x - x_v_eq(1)) ./ (x_v_eq(3) - x_v_eq(1)) .* (x - x_v_eq(2)) ./ (x_v_eq(3) - x_v_eq(2));

figure
% plotto la funzione e i polinomi base
plot(x, f(x), 'LineWidth', 2)
hold on
plot(x, Phi1_eq(x), x, Phi2_eq(x), x, Phi3_eq(x))
plot(x_v_eq, 0*x_v_eq, '*')
plot(x_v_eq, 1, 'k*')

% costruisco il polinomio interpolante di lagrange e plotto
p_eq = @(x) f(x_v_eq(1)) .* Phi1_eq(x) + f(x_v_eq(2)) .* Phi2_eq(x) + f(x_v_eq(3)) .* Phi3_eq(x);

plot(x, p_eq(x), 'LineWidth', 2)

%% verifica

% uso la funzione built-in

P_eq = polyfit(x_v_eq, f(x_v_eq), 2);
P_dis_eq = polyval(P_eq, x);

figure
plot(x, P_dis_eq, 'LineWidth', 2)
hold on
plot(x, p_eq(x))

%% 3

% essendo il polinomio di grado 3 e la funzione da approssimare un
% polinomio di grado 3 troveremo la funzione polinomiale esatta

% verifica
x_v_gr3 = [0, exp(-sqrt(2)), 3^(-sqrt(0.5)), 2];
P_gr3 = polyfit(x_v_gr3, f(x_v_gr3), 3);
P_dis_gr3 = polyval(P_gr3, x);

figure
plot(x, P_dis_gr3, 'LineWidth', 2)
hold on
plot(x, f(x))
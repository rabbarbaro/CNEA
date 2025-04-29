clear
clc
close all

% definisco segnale fisico e rilevazione rumorosa
g = @(x) 10 * x.^2;
f = @(x) g(x) + 2*rand(size(x)) - 1;

% definisco l'intervallo e creo il vettore di query
a = 0;
b = 1;
x_q = linspace(a, b, 1000);

%% 1

% definisco il grado del polinomio interpolante
n = 9;

% definisco n+1 nodi equispaziati e valuto la funzione in quei nodi
x_nod = linspace(a, b, n+1);
f_nod = f(x_nod);
% calcolo il polinomio interpolante di lagrange interpolante sugli n+1 nodi
% definiti prima
P_9 = polyfit(x_nod, f_nod, n);
P_dis_9 = polyval(P_9, x_q);

% definisco il grado del polinomio approssimante ai minimi quadrati
m = 2;

% calcolo il polinomio di grado m approssimante ai minimi quadrati
P_MQ = polyfit(x_nod, f_nod, m);
P_dis_MQ = polyval(P_MQ, x_q);

% plotto il segnale fisico, la rilevazione rumorosa e le due
% approssimazioni trovate
plot(x_q, g(x_q), 'LineWidth', 2)
hold on
plot(x_q, f(x_q))
plot(x_q, P_dis_9, 'LineWidth', 2)
plot(x_q, P_dis_MQ, 'LineWidth', 2)
legend('g(x)', 'f(x)', 'Polinomio gr 9', 'Minimi quadrati gr 2')

%% 2

% estrapolo il segnale in x = 2
x = 2;
valTrue = g(x);
val_P_9 = polyval(P_9, x);
val_P_MQ = polyval(P_MQ, x);

fprintf("Valore estrapolato in x = %f (%f) da:\n" ...
    + "\tPolinomio di Lagrange di grado 9: %f\n" ...
    + "\tPolinomio ai minimi quadrati di grado 2: %f\n", ...
    x, valTrue, val_P_9, val_P_MQ);

% il polinomio di lagrange esplode, non può essere usato per estrapolare
% fuori dall'intervallo in cui sono dati i nodi

%% 3

% l'approssimazione ai minimi quadrati è poco sensibile alle perturbazioni
% casuali, mentre l'interpolazione polinomiale lo è parecchio
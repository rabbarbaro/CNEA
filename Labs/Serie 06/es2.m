clear
clc
close all

% definisco funzione come da testo
f = @(x) sin(1 ./ (1 + x.^2));

% definisco estremi dell'intervallo, lo creo con linspace, valuto la
% funzione sull'intervallo
a = -2*pi;
b = 2*pi;
x = linspace(a, b, 1000);
fx = f(x);

% plotto
plot(x, fx)
grid on

%% 1-2

% creo il vettore con i vari gradi
nv = [2, 4, 8, 10];
% inizializzo le matrici dei polinomi valutati e dell'errore massimo (la
% norma infinito)
PP_dis = [];
err_max = [];

for n = nv
    % creo n+1 nnodi equispaziati
    x_nod = linspace(a, b, n+1);
        % equivalente a fare
        % h = (b - a)/n;
        % x_nod = a:h:b;
    % valuto la funzione nei nodi
    f_nod = f(x_nod);
    % ottengo il polimonio di lagrange dai nodi ottenuti e lo valuto
    % sull'intervallo definito prima 
    P = polyfit(x_nod, f_nod, n);
    P_dis = polyval(P, x);
    % salvo il polinomio valutato e l'errore massimo nelle matrici
    PP_dis = [PP_dis; P_dis];
    err_max = [err_max; max(abs(P_dis - fx))];
end

% plotto funzione e approssimazioni polinomiali
figure
plot(x, fx, x, PP_dis)
grid on
legend('f(x)', 'P2', 'P4', 'P8', 'P10')
title('Approssimazioni polinomiali su nodi equispaziati')

% plotto l'errore rispetto al grado polinomiale
figure
plot(nv, err_max)

% fenomeno di runge!

%% 3

% inizializzo le matrici dei polinomi valutati e dell'errore massimo (norma
% infinito) per i nodi di Chebyshev
PP_dis_che = [];
err_max_che = [];

for n = nv
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
    % salvo il polinomio valutato e l'errore massimo nelle matrici
    PP_dis_che = [PP_dis_che; P_dis_che];
    err_max_che = [err_max_che; max(abs(P_dis_che - fx))];
end

% plotto funzione e approssimazioni polinomiali
figure
plot(x, fx, x, PP_dis_che)      % passando una matrice plotta ogni riga
grid on
legend('f(x)', 'P2', 'P4', 'P8', 'P10')
title('Approssimazioni polinomiali su nodi di Chebyshev')

% plotto l'errore rispetto al grado polinomiale
figure
plot(nv, err_max_che)

% non si presenta fenomeno di runge
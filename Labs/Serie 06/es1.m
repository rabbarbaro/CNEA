clear
clc
close all

%% 1

% definiamo la funzione come da testo
f = @(x) x .* sin(x);

% definiamo gli estremi, creo un vettore equispaziato per plottare
a = -2;
b = 6;
x = linspace(a, b, 1000);
% valuto già la funzione nei vari punti e la salvo in un array
fx = f(x);

% plotto
plot(x, fx)
grid on

%% 2-3

% definisco in un vettore il grado del polinomio approssimante
nv = [2, 4, 6];
% inizializzo le matrici per i polniomi ottenuti (su tutto l'intervallo,
% non solo i coefficienti), l'errore e il massimo dell'errore (la norma
% infinito)
PP_dis = [];
err_dis = [];
err_max = [];

% ciclo for su tutti i gradi
for n = nv
    % creo n+1 nodi equistaziati
    h = (b - a)/n;
    x_nod = a:h:b;
    % calcolo la funzione nei nodi e salvo in un array
    f_nod = f(x_nod);
    % ottengo il polimonio di lagrange dai nodi ottenuti e lo valuto
    % sull'intervallo definito prima
    P = polyfit(x_nod, f_nod, n);
    P_dis = polyval(P, x);
    % salvo il polinomio valutato nella matrice
    PP_dis = [PP_dis; P_dis];
    % calcolo il vettore dell'errore su tutto l'intervallo e lo salvo nella
    % matrice
    err_dis = [err_dis; abs(P_dis - fx)];
    % calcolo il massimo dell'errore e lo salvo nella matrice
    err_max = [err_max; max(abs(P_dis - fx))];
end

figure
plot(x, fx, x, PP_dis)      % passando una matrice plotta per riga
grid on
legend('f(x)', 'P2', 'P4', 'P6')

figure
plot(x, err_dis)
legend('err P2', 'err P4', 'err P6')

idx = 1;
for n = nv
    fprintf("L'errore massimo per il polinomio di grado %d è: %f \n", n, err_max(idx))
    idx = idx + 1;
end
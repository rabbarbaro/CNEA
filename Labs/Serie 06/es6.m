clear
clc
close all

% definisco il vettore delle misure e l'intervallo di tempo
y_s = [10.0 9.89 9.75 9.66 9.10 8.95 8.10 7.49 6.80 6.13 5.05];
deltaT = 0.1;
bigT = 1;
t = 0:deltaT:bigT;

% definisco il vettore del tempo di query
t_q = linspace(0, bigT, 1000);

% definisco la quota iniziale
y0 = 10;

% definisco l'accelerazione di gravità
g = 9.81;

%% 1

% definisco la legge ideale
y = @(t) y0 - 0.5 * g .* t.^2;

% plotto la legge ideale e i dati sperimentali
plot(t_q, y(t_q), t, y_s, '*')

%% 2

% ho il vettore dei nodi (t) e i dati, trovo l'interpolante di lagrange
P = polyfit(t, y_s, length(t)-1);
P_dis = polyval(P, t_q);

% trovo l'interpolante lineare composita
linComp = interp1(t, y_s, t_q);

% cerco l'approssimazione polinomiale ai minimi quadrati di grado 2
P2 = polyfit(t, y_s, 2);
P2_dis = polyval(P2, t_q);

plot(t_q, P_dis, t_q, linComp, t_q, P2_dis, t, y_s, '*')
legend('Polinomio Lagrange', 'Interpolante lineare composito', ...
    'Approssimazione polinomiale di grado 2', 'Dati sperimentali')

% calcolo l'errore
err_P = abs(P_dis - y(t_q));
err_linComp = abs(linComp - y(t_q));
err_P2 = abs(P2_dis - y(t_q));

% plotto l'errore
figure
plot(t_q, err_P, t_q, err_linComp, t_q, err_P2)
legend('Errore Polinomio Lagrange', 'Errore Interpolante lineare composito', ...
    'Errore Approssimazione polinomiale di grado 2')

% trovo il massimo dell'errore (norma infinito)
err_max_P = max(err_P);
err_max_linComp = max(err_linComp);
err_max_P2 = max(err_P2);
fprintf('Errore in norma infinito polinomio Lagrange (massimo): %f\n', err_max_P)
fprintf('Errore in norma infinito interpolante lineare composito (massimo): %f\n', err_max_linComp)
fprintf('Errore in norma infinito approssimazione polinomiale di grado 2 (massimo): %f\n', err_max_P2)

%% 3

% definisco il tempo di query per il calcolo dell'altezza
% e calcolo l'altezza ideale
t_1 = 1.05;
h_ex = y(t_1);
h_P =  polyval(P, t_1);
h_P2 = polyval(P2, t_1);

fprintf('Altezza ideale: %f\n', h_ex)
fprintf('Altezza polinomio Lagrange: %f\n', h_P)
fprintf('Altezza approssimazione polinomiale di grado 2: %f\n', h_P2)

% osserviamo che il polinomio di lagrange esplode
% non è un valido metodo per estrapolare valori fuori dall'intervallo in
% cui sono definiti i nodi
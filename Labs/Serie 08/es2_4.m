clear
clc
close all

%% 1

% scrivo la funzione che implementa il metodo di Crank-Nicolson per la
% risoluzione di un generico problema di Cauchy
% vedi Crank_Nicolson.m

%% 2

% definisco i parametri per il problema di cauchy (problema modello)
y0 = 2;
lambda = -42;

% definisco l'intervallo di tempo
t0 = 0;
t_max = 1;

% definisco la funzione f del problema di cauchy e la soluzione esatta
f = @(t, y) lambda * y;
y_ex = @(t) y0*exp(lambda*(t - t0));

% scelgo il passo e calcolo la soluzione approssimata con Crank-Nicolson
h = 0.02;
[t_h_CN, u_h_CN, iter_pf] = Crank_Nicolson(f, t_max, y0, h);
% plotto soluzione esatta e approssimazione con CN
plot(t_h_CN, y_ex(t_h_CN), t_h_CN, u_h_CN)
legend('Soluzione esatta', 'Approssimazione Crank-Nicolson')

%% 3

% plotto il numero di iterazioni necessarie per il metodo delle iterazioini
% di punto fisso a ogni passo di tempo di CN
figure
plot(t_h_CN, iter_pf, '*-')

%% 4

% scrivo la funzione che implementa il metodo di HEUN per la risoluzione
% di un generico problema di Cauchy
% vedi Heun.m

%% 5

% calcolo la saoluzione approssimata con Heun (h = 0.02)
[t_h_H, u_h_H] = Heun(f, t_max, y0, h);
% plotto soluzione esatta e approssimazione con CH
figure
plot(t_h_H, y_ex(t_h_H), t_h_H, u_h_H)
legend('Soluzione esatta', 'Approssimazione Heun')

%% 6

% definisco il vettore di passi temporali
H = 0.02 * 2.^(0:-1:-4);

% inizializzo i vettori degli errori
e_h_EI = [];
e_h_CN = [];
e_h_H = [];
e_h_RK4 = [];

% per ogni passo e per ogni metodo calcolo l'errore massimo
for h = H
    [t_h_EI, u_h_EI, ~] = eulero_indietro_pto_fisso(f, t_max, y0, h);
    e_h_EI = [e_h_EI, max(abs(u_h_EI - y_ex(t_h_EI)))];

    [t_h_CN, u_h_CN, ~] = Crank_Nicolson(f, t_max, y0, h);
    e_h_CN = [e_h_CN, max(abs(u_h_CN - y_ex(t_h_CN)))];

    [t_h_H, u_h_H] = Heun(f, t_max, y0, h);
    e_h_H = [e_h_H, max(abs(u_h_H - y_ex(t_h_H)))];
    
    [t_h_RK4, u_h_RK4] = Runge_Kutta_4(f, t_max, y0, h);
    e_h_RK4 = [e_h_RK4, max(abs(u_h_RK4 - y_ex(t_h_RK4)))];
end

% plotto in scala logaritmica gli errori trovati
figure
loglog(H, e_h_EI, H, e_h_CN, H, e_h_H, H, e_h_RK4)
hold on
loglog(H, H, H, H.^2, H, H.^4, 'LineStyle','--')
legend('Eulero indietro', 'Crank-Nicolson', 'Heun', 'Runge-Kutta 4', ...
    'h^1', 'h^2', 'h^4')

% come da teoria troviamo:
%   - EI che converge con ordine 1
%   - CN e H che convergono con ordine 2
%   - RK4 che converge con ordine 4
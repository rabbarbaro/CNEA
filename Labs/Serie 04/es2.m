clear
clc
close all

%% 1

% prendiamo la funzione dell'esercizio precedente e la deriviamo (a mano)
f = @(x) x.^3 - (2 + exp(1))*x.^2 + (2*exp(1) + 1)*x + (1 - exp(1)) - cosh(x - 1);
df = @(x) 3*x.^2 - 2*(2 + exp(1))*x + (2*exp(1) + 1) - sinh(x - 1);
% definiamo l'intervallo in cui considerarla
x = 0.5:0.01:6.5;

% plottiamo la funzione, la derivata ed evidenziamo i loro zeri
plot(x, f(x), x, df(x))
grid on
hold on
plot(x, 0*x, '--')
legend("f(x)", "df(x)", '')

% calcolo la derivata seconda della funzione e la valuto in 1
d2f = @(x) 6*x - 2*(2 + exp(1)) - cosh(x - 1);
d2f(1);

% posso quindi applicare il metodo di newton modificato con molteplicictà 
% m = 2 per il primo zero in x = 1 (funzione nulla, derivata prima nulla,
% derivata seconda NON nulla)

%% 2

% creiamo la funzione che implemennta il metodo di newton (modificato)
% vedi file newton.m

%% 3

% definisco la tolleranza come da testo, un numero massimo di iterazioni
% arbitrario e una guess iniziale "abbastanza vicina"
toll = 1e-6;
nmax = 1000;
x0 = .5;

% calcoliamo il primo zero con il metodo normale e quello modificato

[xvect_1, it_1] = newton(x0, nmax, toll, f, df);
fprintf("Il metodo di Newton normale converge allo zero %f in %d iterazioni \n", ...
    xvect_1(end), it_1)

mol = 2;
[xvect_1_m, it_1_m] = newton(x0, nmax, toll, f, df, mol);
fprintf("Il metodo di Newton modificato converge allo zero %f in %d iterazioni \n", ...
    xvect_1_m(end), it_1_m)

% calcoliamo gli altri zeri con newton

x0 = 3;
[xvect_2, it_2] = newton(x0, nmax, toll, f, df);
fprintf("Il metodo di Newton normale converge allo zero %f in %d iterazioni \n", ...
    xvect_2(end), it_2)

x0 = 6;
[xvect_3, it_3] = newton(x0, nmax, toll, f, df);
fprintf("Il metodo di Newton normale converge allo zero %f in %d iterazioni \n", ...
    xvect_3(end), it_3)

% conosciamo il primo zero esatto
alpha1_ex = 1;

% calcoliamo l'errore assoluto per ogni iterazione per il metodo di newton
% normale e per quello modificato sul primo zero

sol_ex_1 = alpha1_ex * ones(size(xvect_1));
err_1 = abs(xvect_1 - sol_ex_1);

sol_ex_1_m = alpha1_ex * ones(size(xvect_1_m));
err_1_m = abs(xvect_1_m - sol_ex_1_m);

% creiamo un vettore da 0 al numero di iterazioni per il metodo di newton
% normale e per quello modificato
it_1v = 0:it_1;
it_1_mv = 0:it_1_m;

% plottiamo l'andamento in scala semilogaritmica dell'errore assoluto per
% il metodo normale e quello modificato
figure
semilogy(it_1v, err_1)
hold on
semilogy(it_1_mv, err_1_m)
grid on
legend("Newton", "Newton modificato")
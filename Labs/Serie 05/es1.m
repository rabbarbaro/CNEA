clear
clc
close all

% definisco la funzione e l'intervallo come da testo
phi = @(y) -y.^4/2 + 4*y.^3 - 7*y.^2 - 4*y + 10;
a = -1;
b = 3;
y = linspace(a, b, 1000);

%% 1

% plotto e verifico esistenza di un minimo in x = 2
x = 2;
plot(y, phi(y), x, phi(x), 'rx')
legend("Funzione obiettivo","Minimo")
xlabel("y")
ylabel("\Phi")

%% 2

% implemento il metodo per approssimare lo zero con la sezione aurea
% vedi sezione_aurea.m

%% 3

% definisco numero massimo di iterazioni e tolleranza
tol = 1e-6;
nmax = 100;

% lancio il metodo della sezione aurea
[xv, k, est_errv] = sezione_aurea(phi, a, b, tol, nmax);

%% 4

% calcolo l'errore effettivo usando come soluzione l'ultima iterazione del
% metodo
errv = abs(x*ones(size(xv)) - xv);

% plotto l'errore e la stima dell'errore
semilogy(est_errv)
hold on
semilogy(errv)
legend('Stima', 'Errore')

% err_est = (b - a)/(2 * ((1 + sqrt(5))/2)^k);
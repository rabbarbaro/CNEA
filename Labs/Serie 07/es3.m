clear
clc
close all
%%
% definisco le funzioni assegnate
p1 = @(x) x.^4 - 2*x + 1;
p2 = @(x) 3*x.^9 - 5*x.^4 + 1;
g = @(x) 10 ./ (x + 2);
z = @(x) sqrt(x);

% definisco gli integrali esatti
Iex_p1 = 1/5;
Iex_p2 = 3/10;
Iex_g = 10 * log(3/2);
Iex_z = 2/3;

% definisco l'intervallo di integrazione
a = 0;
b = 1;

%% 1

% useremo la funzione zplege, che implementa la formula di quadratura
% gaussiana con nodi di Gauss-Legendre

% definisco l'ordine della formula e il numero di nodi
%   - n indica l’ordine della formula
%   - (n + 1 sono il numero di nodi di quadratura)
N = 0:7;

% inizializzo i vettori per gli errori
E_p1 = [];
E_p2 = [];
E_g = [];
E_z = [];

% per ogni ordine e per ogni funzione calcolo l'integrale e l'errore
for n = N
    [x, w] = zplege(n, a, b);

    I_p1 = sum(w .* p1(x));
    I_p2 = sum(w .* p2(x));
    I_g = sum(w .* g(x));
    I_z = sum(w .* z(x));

    E_p1 = [E_p1 abs(I_p1 - Iex_p1)];
    E_p2 = [E_p2 abs(I_p2 - Iex_p2)];
    E_g = [E_g abs(I_g - Iex_g)];
    E_z = [E_z abs(I_z - Iex_z)];
end

%% 2

% plotto gli errori in scala logaritmica
figure
semilogy(N, E_p1, N, E_p2, N, E_g, N, E_z)
legend('p1', 'p2', 'g', 'z')

% troviamo i risultati teorici: per le funzioni polinomiali, ritroviamo che
% le formule di quadratura di GL semplici hanno ordine di convergenza 2n+1
%   - l'errore (analiticamente, poi numericamente fa casino) per la
%     funzione polinomiale p1 di grado 4 è nullo per n = 2 (2n+1 = 5)
%   - l'errore (analiticamente) per la funzione polinomiale p2 di grado 9 è
%     nullo per n = 4 (2n+1 = 9)

%% 3.p1

% calcolo gli integrali con le formule semplici
I_pm = pmedcomp(a, b, 1, p1);
I_tr = trapcomp(a, b, 1, p1);
I_s = simpcomp(a, b, 1, p1);
% calcolo gli errori
E_pm = abs(I_pm - Iex_p1);
E_tr = abs(I_tr - Iex_p1);
E_s = abs(I_s - Iex_p1);

% il confronto con la formula di Gauss-Legendre è legittimo, in quanto
% aumentare il grado della formula di quadratura di Gauss-Legendre non
% aumenta il numero di sottointervalli presi in considerazione
% "esternamente"

% plotto
figure
semilogy(N, E_p1, N, E_pm * ones(size(N)), N, E_tr * ones(size(N)), N, E_s * ones(size(N)))
legend('p1', 'pmedcomp', 'trapcomp', 'simpcomp')

% G-L per n=0 coincide esattamente con punto medio

%% 3.g

% calcolo gli integrali con le formule semplici
I_pm = pmedcomp(a, b, 1, g);
I_tr = trapcomp(a, b, 1, g);
I_s = simpcomp(a, b, 1, g);
% calcolo gli errori
E_pm = abs(I_pm - Iex_g);
E_tr = abs(I_tr - Iex_g);
E_s = abs(I_s - Iex_g);

figure
semilogy(N, E_g, N, E_pm * ones(size(N)), N, E_tr * ones(size(N)), N, E_s * ones(size(N)))
legend('g', 'pmedcomp', 'trapcomp', 'simpcomp')

% G-L per n=0 coincide esattamente con punto medio
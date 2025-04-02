clear
clc
close all

% definisco il funzionale, il suo gradiente e la sua hessiana
Phi = @(y1, y2) 3*y1.^2 + y2.^2 - y1/4 - y2/6 + exp(-2*y1.*y2);
GradPhi = @(y1, y2) [6*y1 - 2*exp(-2*y1.*y2).*y2 - 1/4; ...
                     2*y2 - 2*exp(-2*y1.*y2).*y1 - 1/6];
HessPhi = @(y1, y2) [6 + 4*exp(-2*y1.*y2).*y2^2, exp(-2*y1.*y2).*(4*y1.*y2-2); ...
                     exp(-2*y1.*y2).*(4*y1.*y2-2), 2 + 4*exp(-2*y1.*y2).*y1^2,];

% definisco la soluzione esatta
x_ex = [0.099299729019640; 0.179161952163217];

%% 1

% definisco gli intervalli in y1 e y2 in cui visualizzare la funzione
a1 = -0.15;
b1 = 0.25;
a2 = -0.05;
b2 = 0.35;

% plotto la funzione costo
y1 = linspace(a1, b1, 100);
y2 = linspace(a2, b2, 100);
figure
plot_phi(y1, y2, Phi)

%% 2

% implemento la funzione che esegue il metodo di newton per la ricerca dei
% punti di minimo locali di una funzione costo bidimensionale
% vedi newton_opt.m

%% 3

% stima iniziale
x0 = [-0.14, 0.14]';
tol = 1e-8;
nmax = 200;

% lancio il metodo di newton
[xvect_N, it_N] = newton_opt(GradPhi, HessPhi, x0, tol, nmax);

%% 4

% implemento la funzione che esegue il metodo BFGS per la ricerca dei
% punti di minimo locali di una funzione costo bidimensionale
% vedi bfgs.m

%% 5

% i parametri iniziali sono gli stessi

% lancio il metodo di newton
[xvect_BFGS, it_BFGS] = bfgs(Phi, GradPhi, x0, tol, nmax);

%% 6

% inizializzo i vettori degli errori
err_N = zeros(size(xvect_N, 2), 1);
err_BFGS = zeros(size(xvect_BFGS, 2), 1);

% calcolo vettore errori per newton
for ii = 1:size(xvect_N, 2)
    err_N(ii) = norm(x_ex - xvect_N(:, ii));
end

% calcolo vettore errori per bfgs
for ii = 1:size(xvect_BFGS, 2)
    err_BFGS(ii) = norm(x_ex - xvect_BFGS(:, ii));
end

% plotto
figure
semilogy(err_N,'*-');
hold on;
grid on;
semilogy(err_BFGS,'*-');
legend('Newton','BFGS');
title('Confronto Riduzione errori');
xlabel('Iterazioni');
ylabel('Errore');
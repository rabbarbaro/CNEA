clear
clc
close all

% definisco i coefficienti
C_L0 = 0;
C_Lalpha = 5;
C_Ldelta = 0.3;
C_D0 = 0.02;
C_Dalpha = 0.5;
C_Ddelta = 0.05;

% definisco la funzione e il suo gradiente
Phi = @(y1, y2) -(C_L0 + C_Lalpha*y1 + C_Ldelta*y2) ./ (C_D0 + C_Dalpha*y1.^2 + C_Ddelta*y2.^2);
GradPhi = @(y1, y2) [500*(50*y1.^2 + 6*y1.*y2 - 5*y2.^2 - 2)/(50*y1.^2 + 5*y2.^2 + 2).^2; ...
                     10*(-150*y1.^2 + 500*y1.*y2 + 15*y2.^2 - 6)/(50*y1.^2 + 5*y2.^2 + 2).^2];

% salvo la soluzione esatta
x_ex = [0.196494373584908, 0.117896623783779]';

%% 1

% definisco gli intervalli di y1 e y2 in cui rappresentare la funzione
a1 = 0;
b1 = 0.4;
a2 = -0.1;
b2 = 0.3;

y1 = linspace(a1, b1, 100);
y2 = linspace(a2, b2, 100);

% uso la funzione plot_phi per rappresentare funzione e curve di livello
figure
plot_phi(y1, y2, Phi)

%% 2

% implemento il metodo del gradiente (coniugato) per approssimare un minimo
% locale (metodo di discesa)
% vedi gradiente_coniugato_opt.m

%% 3

% guess iniziale, tolleranza e nunmero massimo di iterazioni
x0 = [0.1, 0.05]';
toll = 1e-5;
nmax = 100;

% calcolo il minimo con il metodo del gradiente e il metodo del gradiente
% coniugato
[xvect_g_1, it_g_1] = gradiente_coniugato_opt(Phi, GradPhi, 'G', x0, toll, nmax);
[xvect_fr_1, it_fr_1] = gradiente_coniugato_opt(Phi, GradPhi, 'FR', x0, toll, nmax);

% plotto soluzioni trovate
figure
plot_phi(y1, y2, Phi, 'Gradiente')
plot_soluzioni(Phi, xvect_g_1, [a1 b1], [a2 b2], x_ex)
figure
plot_phi(y1, y2, Phi, 'Gradiente Coniugato')
plot_soluzioni(Phi, xvect_fr_1, [a1 b1], [a2 b2], x_ex)

%% 4

% cambiamo la guess iniziale
x0 = [0, 0.1]';

% ricalcolo usando metodo del gradiente e del gradiente coniugato
[xvect_g_2, it_g_2] = gradiente_coniugato_opt(Phi, GradPhi, 'G', x0, toll, nmax);
[xvect_fr_2, it_fr_2] = gradiente_coniugato_opt(Phi, GradPhi, 'FR', x0, toll, nmax);

% plotto le soluzioni trovate
figure
plot_phi(y1, y2, Phi, 'Gradiente')
plot_soluzioni(Phi, xvect_g_2, [a1 b1], [a2 b2], x_ex)
figure
plot_phi(y1, y2, Phi, 'Gradiente Coniugato')
plot_soluzioni(Phi, xvect_fr_2, [a1 b1], [a2 b2], x_ex)

% il gradiente coniugato con beta di Fletcher-Reeves non converge in 100
% iterazioni
% cambiando i parametri di backtracking
% c1_bt = 2e-1;
% rho_bt = 0.45;
% nmax_bt = 10;
% effettivamente il gradiente coniugato FR converge in 21 iterazioni

% Verificata la sensibilità dei metodi del gradiente rispetto ai parametri

%% 5

% inizializzo i vettori degli errori
err_g = zeros(size(xvect_g_1, 2), 1);
err_fr = zeros(size(xvect_fr_1, 2), 1);

% calcolo vettore errori per gradiente
for ii = 1:size(xvect_g_1, 2)
    err_g(ii) = norm(x_ex - xvect_g_1(:, ii));
end

% calcolo vettore errori per gradiente coniugato
for ii = 1:size(xvect_fr_1, 2)
    err_fr(ii) = norm(x_ex - xvect_fr_1(:, ii));
end

% plotto
figure
semilogy(err_g,'*-');
hold on;
grid on;
semilogy(err_fr,'*-');
legend('Gradiente','Gradiente Coniugato (FR)');
title('Confronto Riduzione errori');
xlabel('Iterazioni');
ylabel('Errore');
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

% definisco gli intervalli di y1 e y2 in cui rappresentare la funzione
a1 = 0;
b1 = 0.4;
a2 = -0.1;
b2 = 0.3;

y1 = linspace(a1, b1, 100);
y2 = linspace(a2, b2, 100);

%% 1

% implemento il metodo del gradiente (coniugato) per approssimare un minimo
% locale (metodo di discesa) con scelta del parametro alpha con
% backtracking o sezione aurea
% vedi gradiente_coniugato_opt_BT_SA.m

%% 2

% guess iniziale, tolleranza e nunmero massimo di iterazioni
x0 = [0.1, 0.05]';
toll = 1e-5;
nmax = 100;

% calcolo il minimo con il metodo del gradiente e il metodo del gradiente
% coniugato usando la sezione aurea per alpha
[xvect_G_SA, it_G_SA] = gradiente_coniugato_opt_BT_SA(Phi, GradPhi, 'G', 'SA', x0, toll, nmax);
[xvect_FR_SA, it_FR_SA] = gradiente_coniugato_opt_BT_SA(Phi, GradPhi, 'FR', 'SA', x0, toll, nmax);

% Visualizzo le soluzioni
figure()
plot_phi(y1, y2, Phi, 'Gradiente con SA');
plot_soluzioni(Phi, xvect_G_SA, [a1 b1], [a2 b2], x_ex);
figure()
plot_phi(y1, y2, Phi, 'Gradiente coniugato FR con SA');
plot_soluzioni(Phi, xvect_FR_SA, [a1 b1], [a2 b2], x_ex);

%% 3

% calcolo i vettori con gli errori
err_G_SA = zeros(size(xvect_G_SA,2),1);
err_FR_SA = zeros(size(xvect_FR_SA,2),1);
for ii = 1:size(xvect_G_SA,2)
    err_G_SA(ii) = norm(x_ex - xvect_G_SA(:,ii));
end
for ii = 1:size(xvect_FR_SA,2)
    err_FR_SA(ii) = norm(x_ex - xvect_FR_SA(:,ii));
end

% plotto in scala semilogaritmica
figure()
semilogy(0:it_G_SA, err_G_SA,'*-')
hold on;
grid on;
semilogy(0:it_FR_SA, err_FR_SA,'*-');
legend('GradSA','ConjGradFRSA');
title('Confronto riduzione errori');
xlabel('#iterazioni');
ylabel('Errori');
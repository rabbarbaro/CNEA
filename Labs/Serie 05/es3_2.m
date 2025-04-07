clear
clc
close all

% definisco funzione, gradiente ed hessiana
Phi = @(y1, y2) 3*y1.^2 + y2.^2 - y1/4 - y2/6 + exp(-2*y1.*y2);
GradPhi = @(y1, y2) [6*y1 - 2*exp(-2*y1.*y2).*y2 - 1/4; ...
                     2*y2 - 2*exp(-2*y1.*y2).*y1 - 1/6];
HessPhi = @(y1, y2) [6 + 4*exp(-2*y1.*y2).*y2^2, exp(-2*y1.*y2).*(4*y1.*y2-2); ...
                     exp(-2*y1.*y2).*(4*y1.*y2-2), 2 + 4*exp(-2*y1.*y2).*y1^2,];

% definisco la soluzione esatta
x_ex = [0.099299729019640; 0.179161952163217];

% guess iniziale, tolleranza e numero massimo di iterazioni assegnate
x0 = [-0.14, 0.14]';
tol = 1e-8;
nmax = 200;

% risolvo con i 4 metodi
[xvect_N, it_N] = newton_opt(GradPhi, HessPhi, x0, tol, nmax);
[xvect_BFGS, it_BFGS] = bfgs(Phi, GradPhi, x0, tol, nmax);
[xvect_G, it_G] = gradiente_coniugato_opt(Phi, GradPhi, 'G', x0, tol, nmax);
[xvect_FR, it_FR] = es3_2_gradiente_coniugato_opt(Phi, GradPhi, 'FR', x0, tol, nmax);

% calcolo i vettori degli errori
errvec_N = [];
for k = 1:it_N
    errvec_N = [errvec_N, norm(xvect_N(:, k) - x_ex)];
end

errvec_BFGS = [];
for k = 1:it_BFGS
    errvec_BFGS = [errvec_BFGS, norm(xvect_BFGS(:, k) - x_ex)];
end

errvec_G = [];
for k = 1:it_G
    errvec_G = [errvec_G, norm(xvect_G(:, k) - x_ex)];
end

errvec_FR = [];
for k = 1:it_FR
    errvec_FR = [errvec_FR, norm(xvect_FR(:, k) - x_ex)];
end

% calcolo i mu per i vari ordini di convergenza
conv_ord1_N = errvec_N(2:end) ./ errvec_N(1:end-1);
conv_ord2_N = errvec_N(2:end) ./ errvec_N(1:end-1).^2;

conv_ord1_BFGS = errvec_BFGS(2:end) ./ errvec_BFGS(1:end-1);
conv_ord2_BFGS = errvec_BFGS(2:end) ./ errvec_BFGS(1:end-1).^2;

conv_ord1_G = errvec_G(2:end) ./ errvec_G(1:end-1);
conv_ord2_G = errvec_G(2:end) ./ errvec_G(1:end-1).^2;

conv_ord1_FR = errvec_FR(2:end) ./ errvec_FR(1:end-1);
conv_ord2_FR = errvec_FR(2:end) ./ errvec_FR(1:end-1).^2;

% plotto

figure

subplot(2,2,1)
plot(conv_ord1_N,'o-');
hold on;
plot(conv_ord2_N,'o-');
legend('ordine1','ordine2');
xlabel('#iterazioni')
ylabel('\mu')
ylim([0,2]);
title('Newton');


subplot(2,2,2)
plot(conv_ord1_BFGS,'o-');
hold on;
plot(conv_ord2_BFGS,'o-');
legend('ordine1','ordine2');
xlabel('#iterazioni')
ylabel('\mu')
ylim([0,2]);
title('BFGS');


subplot(2,2,3)
plot(conv_ord1_G,'o-');
hold on;
plot(conv_ord2_G,'o-');
legend('ordine1','ordine2');
xlabel('#iterazioni')
ylabel('\mu')
ylim([0,2]);
title('Gradiente');

subplot(2,2,4)
plot(conv_ord1_FR,'o-');
hold on;
plot(conv_ord2_FR,'o-');
legend('ordine1','ordine2');
xlabel('#iterazioni')
ylabel('\mu')
ylim([0,2]);
title('Gradiente coniugato');

plot(conv_ord1_FR,'o-');
hold on;
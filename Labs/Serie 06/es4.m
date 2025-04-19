clear
clc
close all

% salviamo i dati forniti
sigma = [0.1800 0.3000 0.5000 0.6000 0.7200 0.7500 0.8000 0.9000 1.0000];       % [1000 kgf/cm^2]
epsilon = [0.0005 0.0010 0.0013 0.0015 0.0020 0.0045 0.0060 0.0070 0.0085];     % [cm/cm]

% plot dei dati
plot(sigma, epsilon, '*')
grid on

a = sigma(1);
b = sigma(end);
sigma_q = linspace(a, b, 1000);

%% 1.1 (interpolazione polinomiale di lagrange)

% ottengo il grado del polinomio
n = length(sigma) - 1;
% valuto il polinomio approssimante di lagrange
P = polyfit(sigma, epsilon, n);
epsilon_P = polyval(P, sigma_q);

figure
plot(sigma, epsilon, '*')
grid on
hold on
plot(sigma_q, epsilon_P)
legend('Dati', 'Interpolazione polinomiale')

% ESPLODE! non va bene

%% 1.2 (interpolazione polinomiale lineare a tratti)

% calcolo l'interpolazione lineare a tratti
epsilon_ICL = interp1(sigma, epsilon, sigma_q);

figure
plot(sigma, epsilon, '*')
grid on
hold on
plot(sigma_q, epsilon_ICL)
legend('Dati', 'Interpolazione lineare a tratti')

%% 1.3 (splines)

% calcolo la spline cubica naturale
epsilon_S_N = cubicspline(sigma, epsilon, sigma_q);

% calcolo la spline not-a-knot
epsilon_S_NAK = spline(sigma, epsilon, sigma_q);

figure
plot(sigma, epsilon, '*')
grid on
hold on
plot(sigma_q, epsilon_S_N, sigma_q, epsilon_S_NAK)
legend('Dati', 'Spline cubica naturale', 'Spline not-a-knot')

% molto simili, cambiano condizioni di chiusura
% c'è un dip fra 0.6 e 0.7, non è un comportamento accettabile fisicamente
% non va bene

%% 1.4 (approssimazione ai minimi quadrati)

% calcolo i polinomi ai minimi quadrati per i gradi richiesti
p_mq1 = polyfit(sigma, epsilon, 1);
p_mq2 = polyfit(sigma, epsilon, 2);
p_mq4 = polyfit(sigma, epsilon, 4);

% valuto su xq
p_mq1_dis = polyval(p_mq1, sigma_q);
p_mq2_dis = polyval(p_mq2, sigma_q);
p_mq4_dis = polyval(p_mq4, sigma_q);

figure
plot(sigma, epsilon, '*')
grid on
hold on
plot(sigma_q, p_mq1_dis, sigma_q, p_mq2_dis, sigma_q, p_mq4_dis)
legend('Dati', 'MQ1', 'MQ2', 'MQ4')

%% 2

figure
plot(sigma, epsilon, '*')
grid on
hold on
plot(sigma_q, epsilon_P, sigma_q, epsilon_ICL, sigma_q, epsilon_S_N, ...
    sigma_q, epsilon_S_NAK, sigma_q, p_mq4_dis)
ylim([0, 0.01])
legend('Dati', 'Interpolazione polinomiale', 'Interpolazione lineare a tratti', ...
    'Spline cubica naturale', 'Spline not-a-knot', 'MQ4', Location='northwest')


%% 3

% punti di sigma in cui calcolare epsilon
sigmav = [0.4, 0.65];                           % [1000 kgf/cm^2]

% interpolazione polinomiale
epsilon_Pv = polyval(P, sigmav)
% questa non ha un minimo di senso, esplode (valori negativi non hanno
% alcun senso fisico)

% interpolazione lineare a tratti
epsilon_ICLv = interp1(sigma, epsilon, sigmav)
% ha senso

% spline cubica naturale
epsilon_S_Nv = cubicspline(sigma, epsilon, sigmav)
% per 0.65 c'è un dip, non ha senso fisico

% minimi quadrati di grado 4
p_mq4_disv = polyval(p_mq4, sigmav)
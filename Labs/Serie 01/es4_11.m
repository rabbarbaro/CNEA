clear
clc
close all

%% circonferenza

% tentativo scacato
% x = -1:0.01:1;
% y = sqrt(1 - x.^2);
% plot(x, y, 'k', x, -y, 'k');

% in coordinate polari
theta = 0:0.01:2*pi;
rho = 1;
x = rho * cos(theta);
y = rho * sin(theta);
plot(x,y)
axis equal
grid on
hold on

%% ellisse

% anche l'ellisse in polari
a = 2;
b = 3;
x = a * cos(theta);
y = b * sin(theta);
plot(x,y)

%% quadrato

% per il quadrato definisco 5 punti (primo e ultimo coincidono)
plot([-2.5 -2.5 2.5 2.5 -2.5], [-2.5, 2.5, 2.5, -2.5, -2.5])

legend('Circonferenza', 'Ellisse', 'Quadrato')
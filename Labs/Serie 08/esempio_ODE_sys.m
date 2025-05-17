clear
clc
close all

% definisco istante iniziale, finale e condizioni inziali
t0 = 0;
tf = 10;
y0 = [2, 2];

% chiamo la risoluzione del sistema di ODE fun_problema
[t, u] = ode45(@fun_problema, [t0, tf], y0);

% plotto la soluzione (entrambe le variabili)
plot(t, u(:, 1), t, u(:, 2))

%% ODE functions

% funzione che definisce il sistema di ODE
function fn = fun_problema(t, y)
    [n, m] = size(y);
    fn = zeros(n, m);       % non strettamente necessario (?)
    % definisco il sistema di ODE
    fn(1) = y(1) - (1 + sin(pi * t)) * y(1) * y(2);
    fn(2) = y(2) * y(1) - y(2);
end
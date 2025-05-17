clear
clc
close all

%% 1

% vedi funzione mms in fondo al file

%% 2

% definisco istante iniziale, finale e condizioni inziali
t0 = 0;
tf = 100;
y0 = [2, 0];

% chiamo la risoluzione del sistema di ODE mms
[t, u] = ode45(@mms, [t0, tf], y0);

% definisco la funzione di un'oscillazione smorzata con i parametri
% corretti
gamma = 0.1;
omega = 1;
omega_1 = sqrt(4*omega^2 - gamma^2)/2;
y_dmp_osc = @(t) 2*exp(-0.05 * t) .* cos(omega_1 * t);

% plotto la soluzione (posizione x = u1), confrontando con oscillazione
% smorzata
plot(t, u(:, 1), t, y_dmp_osc(t), '--')
legend('Soluzione numerica', 'Soluzione esatta')

%% 3

% chiamo la risoluzione del sistema di ODE mms_forz
[t_f, u_f] = ode45(@mms_forz, [t0, tf], y0);

% definisco la funzione di un'oscillazione armonica
omegaf = 0.5;
y_har = @(t) 0.665 * sin(omegaf * t);

plot(t_f, u_f(:, 1), t_f, y_har(t_f), '--')
legend('Soluzione numerica', 'Soluzione a regime')

%% ODE functions

% funzione che definisce il sistema di ODE (non forzato)
function fn = mms(t, y)
    gamma = 0.1;
    omega = 1;

    % x'' + gamma*x' + omega^2*x = 0
    % x'' = - gamma*x' - omega^2*x

    % y(1) = x
    % y(2) = x'

    % y'(1) = y(2)
    % y'(2) = - gamma*y(2) - omega^2*y(1)

    [n, m] = size(y);
    fn = zeros(n, m);
    % definisco il sistema di ODE
    fn(1) = y(2);
    fn(2) = - gamma*y(2) - omega^2*y(1);
end

% funzione che definisce il sistema di ODE (forzato)
function fn = mms_forz(t, y)
    gamma = 0.1;
    omega = 1;
    A0 = 0.5;
    omegaf = 0.5;

    % x'' + gamma*x' + omega^2*x = A0*sin(omegaf * t)
    % x'' = - gamma*x' - omega^2*x + A0*sin(omegaf * t)

    % y(1) = x
    % y(2) = x'

    % y'(1) = y(2)
    % y'(2) = - gamma*y(2) - omega^2*y(1) + A0*sin(omegaf * t)

    [n, m] = size(y);
    fn = zeros(n, m);
    % definisco il sistema di ODE
    fn(1) = y(2);
    fn(2) = - gamma*y(2) - omega^2*y(1) + A0*sin(omegaf * t);
end
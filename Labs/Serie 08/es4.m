clear
clc
close all

% definisco istante iniziale e finale
t0 = 0;
tf = 10;

% dato il problema di Cauchy con ODE del secondo ordine:
%       x'' + 2*x' + 10*x^2 = z
%       x'_0 = 2
%       x_0 = 0
% dove:
%       z = exp(-t/2) .* (2*cos(t) -7/2 * sin(t)) + 40*exp(-t).*(sin(t))^2

%% 1

% il problema è riscrivibile come sistema di ODE del primo ordine:
%       x' = w
%       w' = - 2*w - 10*x^2 + z
%       w_0 = 2
%       x_0 = 0
% nella notazione compatta di sistema:
%       y_v' = f_v + g_v
%       y_v(0) = y_v_0
% in cui il vettore y_v è definito come:
%       y_v = [x, w]'

% scrivendo i vettori g, f, y0: 

g = @(t) [0;
          exp(-t/2) .* (2*cos(t) -7/2 * sin(t)) + 40*exp(-t).*(sin(t))^2];

f = @(y) [y(2);
          -2*y(2) - 10*y(1).^2];

y0 = [0, 2]';

%% 2

% definisco passo h = 0.01
h = 1e-2;
% calcolo il numero di passi temporali (h = (tf-t0)/Nh)
Nh = round((tf-t0)/h);
% genero l'intervallo temporale dagli estremi
tv = [t0, tf];

% salvo come funzione vettoriale in t e y il right-hand-side del sistema
fun = @(t, y) f(y) + g(t);

% chiamo la funzione eulero_avanti_sistemi
[t, u] = eulero_avanti_sistemi(fun, tv, y0, Nh);

% plotto la soluzione (x e x')
plot(t, u(1, :), t, u(2, :))
legend('x', "x'")

% valuto l'approssimazione di x(t_1), x(t_f), con t_n = nh, h = 0,1,...,Nh
u(1, 2)         % dato che il primo elemento del vettore è la CI
u(1, end)

% metodo grafico = zoom su grafico --> data tips per scorrere indici
plot(t, abs(u(1, :)), t, 0.1*ones(1, length(t)))
legend('x', "0.1")
% trovo 5.68

% alternativamente, metodo automatico con find
idx = find((abs(u(1, :)) - 0.1) > 0, 1, 'last');
t(idx+1)

% terzo metodo, dove cerco iterativamente per ogni valore di k (indice sul
% tempo) se il modulo della soluzione scende sotto a 0.1, e poi prendo
% l'ultimo di questi valori trovati
tm = [];
for k = 2:length(t)
    if abs(u(1, k-1)) >= 0.1 && abs(u(1, k)) < 0.1
        tm = [tm, t(k)];
    end
end
tm = tm(end)

%% 3

% definisco la soluzione esatta data dal testo
x_ex = @(t) 2*exp(-t/2).*sin(t);

% definisco il vettore dei passi
hv = 1e-3 ./ 2.^(0:3);

% inizializzo il vettore degli errori
err_EA = [];

% per ogni passo
for h = hv
    % calcolo il numero di istanti discreti
    Nh = round((tf-t0)/h);
    % chiamo la soluzione con il metodo di eulero in avanti
    [t, u] = eulero_avanti_sistemi(fun, tv, y0, Nh);
    % calcolo l'errore massimo per ogni passo
    err_EA = [err_EA, max(abs(u(1,:) - x_ex(t)))];
end

%% 4

% sapendo dalla teoria che:
%       E <= C*h^p
%       Eh1 / Eh2 = (h1 / h2)^p
%       p = log(Eh1 / Eh2) / log(h1 / h2)

p = log(err_EA(2:end)./err_EA(1:end-1)) ./ log(hv(2:end) ./ hv(1:end-1))

%% 5

% conosco la soluzione esatta, quindi posso calcolare la matrice jacobiana
% della funzione f rispetto alla variabile y per ogni passo
% J(t) = [0, 1;
%         -20*y(1), -2];
% J(t) = [0, 1;
%         -40*exp(-t/2).*sin(t), -2]

J = @(t) [0, 1;
          -20*x_ex(t), -2];

% per avere asintotica stabilità, ogni autovalore di J deve avere parte
% reale negativa
% per avere assoluta stabilità del metodo numerico devo garantire che ogni
% autovalore resti dentro alla regione di assoluta stabilità del metodo
% nel caso di EA la formula chiusa è:
%       hmax = - 2 * (real(lambda(i) / abs(lambda(i)^2)))

% creo il vettore di tempo su cui calcolare il massim degli autovalori
tv = t0:0.01:tf;
% inizializzo i vettori per gli autovalori e per l'h calcolato
lambdav = [];
hv = [];

% per ogni istante di tempo:
for t = tv
    % calcolo gli autovalori e li salvo nel vettore
    lambda = eig(J(t));
    lambdav = [lambdav, lambda];
    % calcolo l'h max corrispondente per EA e lo salvo nel vettore
    hmax = - 2 * (real(lambda(2) / abs(lambda(2)^2)));
    hv = [hv, hmax];
end

% calcolo il minimo degli hmax
hmax = min(hv)

% per visualizzare:
%   - il parametro h è un riscalamento dell'autovalore visto come vettore
%   - devo fare in modo che l'autovalore con modulo maggiore (e parte reale
%     negativa) resti dentro alla regione di assoluta stabilità del metodo
h = 1;
plot(real(h*lambdav), imag(h*lambdav), 'xr')
hold on
plot(real(hmax*lambdav), imag(hmax*lambdav), 'xb')
xline(0, 'k--');
theta = linspace(0, 2*pi, 200);
plot(-1 + cos(theta), sin(theta), 'k-')
legend('Autovalori (lambda * h con h = 1)', 'lambda * hmax', 'Re = 0', 'Regione di assoluta stabilità EA')
hold off
grid on
axis equal

% ricordando che la regione di stabilità del metodo di eulero in avanti
% è il cerchio unitario centrato in -1
% e che ha senso cercare l'assoluta stabilità del metodo numerico solo se
% il problema matematico accetta una soluzione asintoticamente stabile,
% ovvero solo se gli autovalori hanno parte reale negativa

% nel caso non avessi la forma chiusa:
%       h * lambda_max = (h * real(lambda_max), h * imag(lambda_max))
% dove lambda_max è l'autovalore più lontano da -1, poi imponendo:
%       |1 + h * lambda_max| < 1
% e poi risolvendo la disequazione risultante

% in generale vale |R(z)| < 1 per qualsiasi metodo, con R la funzione di
% stabilià

%% 6

% implemento il metodo multipasso come da testo
function [t, u] = multipasso(fun, tv, y0, Nh)

h = (tv(end) - tv(1)) / Nh;
t = linspace(tv(1), tv(end), Nh+1);
u = zeros(size(y0, 1), Nh);

u(:, 1) = y0;

u(:, 2) = u(:, 1) + h * fun(t(1), u(:, 1));

for n = 2:Nh
    u(:, n+1) = u(:, n) + 3/2 * h * fun(t(n), u(:,n)) - 1/2 * h * fun(t(n-1), u(:,n-1));
end
end

% definisco passo h = 0.01
h = 1e-2;
% calcolo il numero di passi temporali (h = (tf-t0)/Nh)
Nh = round((tf-t0)/h);
% genero l'intervallo temporale dagli estremi
tv = [t0, tf];

% chiamo il metodo multipasso
[t, u] = multipasso(fun, tv, y0, Nh);

% plotto la soluzione (x e x')
plot(t, u(1, :), t, u(2, :))
legend('x', "x'")

% valuto l'approssimazione di x(t_1), x(t_1), x(t_f), con t_n = nh,
% h = 0,1,2...,Nh
u(1, 2)         % dato che il primo elemento del vettore è la CI
u(1, 3)
u(1, end)

%% 7

% definisco il vettore dei passi
hv = 1e-3 ./ 2.^(0:3);
% inizializzo il vettore degli errori
err_MS = [];

% per ogni passo
for h = hv
    % calcolo il numero di istanti discreti
    Nh = round((tf-t0)/h);
    % chiamo la soluzione con il metodo di eulero in avanti
    [t, u] = multipasso(fun, tv, y0, Nh);
    % calcolo l'errore massimo per ogni passo
    err_MS = [err_MS, max(abs(u(1,:) - x_ex(t)))];
end

% sapendo dalla teoria che:
%       E <= C*h^p
%       Eh1 / Eh2 = (h1 / h2)^p
%       p = log(Eh1 / Eh2) / log(h1 / h2)

p = log(err_MS(2:end)./err_MS(1:end-1)) ./ log(hv(2:end) ./ hv(1:end-1))
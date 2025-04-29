clear
clc
close all

% minimi quadrati non lineari

% campo di velocità di un getto turbolento libero mediato in tempo
% (Reynolds-averaged)
% la componente radiale è trascurabile, può essere descritto attraverso
% la sola componente assiale come:
%    u(a,r) = (U*d)/a * exp(-A*(r/a)^2)
% con:
%   - a, distanza assiale dal punto di uscita del getto (ugello)
%   - r, distanza radiale dal punto di uscita del getto (ugello)
%   - U, velocità caratteristica del getto
%   - A, parametro adimensionale che determina la decrescita in direzione
%       radiale
%   - d, diametro dell'ugello

% abbiamo:
%   - n, numero di coppie di dati
%   - m, numero di parametri da stimare

% volendo minimizzare la somma dei quadrati delle differenze tra i dati e
% il modello, si ha:
%   sum_i=0^n (y_i - f_tilde(x_i; w_vec))^2

% min_(x_vec in R^m) Phi(x_vec)
% Phi(y_vec) = 1/2 ||F_vec(y_vec)||^2 = 1/2 * sum_i=0^n (f_i(x_i))^2
% con il gradiente e l'hessiano definiti come:
%   grad_Phi(y_vec) = J(y_vec)^T * F(y_vec)
%   Hess_Phi(y_vec) = J(y_vec)^T * J(y_vec) + Q(y_vec)

% usiamo questo gradiente e hessiano per calcolare il minimo usando
% l'hessiana semplificata eliminando Q abbiamo il metodo di Gauss-Newton

%% 1

% implementiamo la funzione gauss_newton che calcola il minimo di una
% funzione non lineare usando il metodo di Gauss-Newton
% vedi gauss_newton.m

%% 2

% definisco i risultati veri
U = 1;
A = 20;
% definisco i dati
d = 1;
a = 0.8;

% definisco il vettore di funzioni f (sui vari dati)
f = @(x, w) (d * w(1))/a * exp(-w(2) * (x/a).^2);
% definisco il vettore dei residui (sui vari dati)
r = @(x, y, w) y - f(x, w);
% definisco lo jacobiano (derivo i residui, è una matrice)
J = @(x, y, w) [- (d/a) * exp(-w(2) * (x/a).^2), ...
                (d * w(1))/a * exp(-w(2) * (x/a).^2) .* (x/a).^2];

% calcolo i 50 nodi su cui calcoliamo i dati (deve essere un vettore
% colonna perché f (e quindi r) è un vettore 50*1
x = linspace(-0.5, 0.5, 50)';

% definisco il vettore w con i risultati veri
w_true = [U, A]';

% calcoliamo il valore della funzione con i dati veri per avere il valore
% nei nodi
y = f(x, w_true);

% definisco la guess iniziale, la tolleranza e il massimo numero di
% iterazioni
w0 = [2, 2]';
toll = 1e-6;
maxit = 1000;

% chiamiamo la funzione
[w_vect, it] = gauss_newton(x, y, r, J, w0, toll, maxit);
w_vect(:, end)

% plottiamo
plot(x, y, 'o', x, f(x, w_vect(:, end)));
% notiamo che abbiamo effettivamente ritrovato i coefficienti corretti,
% dato che la legge trovata si sovrappone ai dati

%% 3

% riapplichiamo il metodo di Gauss-Newton con una funzione a cui abbiamo
% aggiunto del noise gaussiano
y_noise = y + 3e-2 * randn(size(x));
[w_vect_noise, it_noise] = gauss_newton(x, y_noise, r, J, w0, toll, maxit);
w_vect_noise(:, end)

% plottiamo
plot(x, y_noise, 'o', x, f(x, w_vect_noise(:, end)));
% i risultati non sono perfetti ma offrono un'approssimazione robusta
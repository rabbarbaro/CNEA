clear
clc
close all

% definisco la funzione vettoriale e la sua jacobiana
F = @(x1, x2) [6*x1 - 2*exp(-2*x1.*x2).*x2 - 1/4; ...
               2*x2 - 2*exp(-2*x1.*x2).*x1 - 1/6];
JF = @(x1, x2) [6 + 4*exp(-2*x1.*x2).*x2^2, exp(-2*x1.*x2).*(4*x1.*x2-2); ...
               exp(-2*x1.*x2).*(4*x1.*x2-2), 2 + 4*exp(-2*x1.*x2).*x1^2,];

% prendo dal testo la soluzione esatta
alpha = [0.099299729019640, 0.179161952163217]';

%% 1

% iterata iniziale
x = [-0.14, 0.14]';
% numero massimo di iterazioni
kmax = 5;

% inizializzo xvec come vettore e l'iterata 1
xvec = [];

% per ogni iterazione
for k = 1:kmax
    % valuto funzione e jacobiana nell'iterata corrente
    Fvec = F(x(1), x(2));
    Jmat = JF(x(1), x(2));
    % risolvo il sistema lineare
    x = x - Jmat\Fvec;
    % salvo l'iterata
    xvec = [xvec x];
end

% effettivamente è vicina alla soluzione esatta
xvec(:, end);

%% 2

% verifico che il determinante della jacobiana sia != 0
det(JF(alpha(1), alpha(2)));
% ci aspettiamo convergenza quadratica

% per stimare l'ordine di convergenza

% calcoliamo il vettore degli errori
errvec = [];
for k = 1:kmax-1
    errvec = [errvec, norm(xvec(:, k) - alpha)];
end

% calcoliamo un vettore contente il rapporto fra errori successivi
% elevando il denominatore all'ordine di convergenza troveremo il fattore
% di convergenza
% aumentando l'ordine di convergenza, il primo che avrà un fattore di
% convergenza non nullo o non infinito sarà l'ordine di convergenza
% dalla teoria:
% lim ||x^(k+1)-alpha||/||x^(k)-alpha||^p = mu

% p = 1?
conv_ord1 = errvec(2:end) ./ errvec(1:end-1)
% la convergenza di ordine 1 tende a 0 --> ordine superiore

% p = 2?
conv_ord2 = errvec(2:end) ./ errvec(1:end-1).^2
% non tende a zero: ordine 2
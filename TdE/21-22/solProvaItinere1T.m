clear
clc
close all

%% 1

beta = 2;
t = 4;
m = bin2dec("1011");
s = 0;
e = 2;

x = (-1)^s * m * beta^(e - t)

%% 2

n = 90;
S = 0;

for k = 0:n
    S = (S + (-1)^k / (2*k + 1));
end
S*4

%% 3

% per ogni riga completa faremo due sottrazioni e due somme
% la prima riga è semplice assegnazione
% la seconda riga è una sottrazione e una moltiplicazione

n = 2000;
nflops = 2 + (n-2)*4

%% 4

n = 7;
A = hilb(n);
b = 5*ones(n, 1);
x = A\b;

res = b - A * x;
err_rel_est = cond(A) * norm(res) / norm(b)

%% 5

% 3.7 + 3/10 * gamma

% TEST
% gamma = 4;
% A = [3 3
%     gamma 1];
% x0 = [0; 1];
% 
% [lambda,x,iter] = eigpower(A,1e-12,1,x0);
% lambda

%% 6

n = 100;
A = 2*diag(ones(n, 1)) - diag(ones(n-1, 1), -1) - diag(ones(n-1, 1), +1);

jj = 47;
lambda_47 = 2 + 2*cos(pi * jj/(n+1));

jj = 46;
lambda_46 = 2 + 2*cos(pi * jj/(n+1));
jj = 48;
lambda_48 = 2 + 2*cos(pi * jj/(n+1));

upper_bound_s = (lambda_47 + lambda_46)/2
lower_bound_s = (lambda_47 + lambda_48)/2
% s != lambda_47

%% 7

% ???

% copilot risponde, con fonte
% https://math.libretexts.org/Bookshelves/Calculus/CLP-1_Differential_Calculus_(Feldman_Rechnitzer_and_Yeager)/06%3A_Appendix/6.03%3A_C-_Root_Finding/6.3.02%3A_C.2_The_Error_Behaviour_of_Newton's_Method

% Definizione della funzione e delle sue derivate
f = @(x) exp(3*x) - 1;
df = @(x) 3*exp(3*x);
d2f = @(x) 9*exp(3*x);

% Valore dell'errore al passo k
error_k = 1e-2;

% Calcolo dell'errore al passo k+1
% Valutiamo f' e f'' nel punto x^{(k)}
xk = 0.1; % Supponiamo xk per il calcolo (valore vicino allo zero della funzione)
f_prime = df(xk);
f_double_prime = d2f(xk);

% Calcolo dell'errore al passo k+1
error_k_plus_1 = (f_double_prime / (2 * f_prime)) * error_k^2;

% Stampa dei risultati
fprintf('Errore al passo k: %.10e\n', error_k);
fprintf('Errore stimato al passo k+1: %.10e\n', error_k_plus_1);

%% 8

f = @(x) exp(-4*x) - 2*x;
x0 = 1;
tol = 1e-2;

err = tol + 1;
k = 0;
while err > tol
    x = x0 - (f(x0)^2) / (f(x0 + f(x0)) - f(x0));
    x0 = x;
    err = abs(f(x));
    k = k+1;
end

k
x

%% 9

phi = @(x) x - 9/2 * log(x / 3);
x0 = 2;

[succ, it] = ptofis(x0, phi, 4, 1e-12);
succ(end)

%% 10

phi = @(x) x - 140/11 * (exp(x/7 - 1) - 1);
dphi = @(x) 1 - 20/11 * exp(x/7 - 1);
alpha = 7;

x = 4:.001:8;

plot(x, phi(x), x, 8*ones(length(x), 1))
axis equal

% |phi'(x)| < 1
% x < -7 * (-1 + log(2) + log(5) - log(11))
% x < 7.6672

% phi'(x) < 0 sempre in [a, b]
phi(4)      % fuori da [a, b]
phi(8)      % dentro [a, b]
% x - 140/11 * (exp(x/7 - 1) - 1) = 8
% per via grafica trovo x = 5.4013
% phi(x) non appartiene ad [a, b] per x < alpha)

% 5.4013 <= x < 7.6672
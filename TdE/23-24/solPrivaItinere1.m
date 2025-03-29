clear
clc
close all

%% T1

% F(2, t -6, 6)
% errore assoluto |x - fl(x)| <= 0.5 * eps_M * x
% eps_M = beta^(1-t)
% 0.5 * 2^(1-t) * x = 2^t * x < err_max
% -t < log_2(err_max / x)

x = sqrt(pi);
err_max = 1e-10;
t = - log2(err_max / x);
t = ceil(t);

%% T2

% A bidiagonale superiore, posso risolvere direttamente i sistemi senza
% fattorizzazione
% 3 * (n - 1) + 1 flops

n = 20;
nSys = 50;
nflops = 50 * (3 * (n - 1) + 1);

%% T3

% la matrice è triangolare superiore --> sostituzioni all'indietro
% err_rel <= K_2(A) * res_rel

n = 10;
res_rel = 1e-12;
A = zeros(n);
for ii = 1:n
    for jj = ii:n
        A(ii, jj) = 21 - ii - jj;
    end
end
err_rel_est = cond(A) * res_rel;

%% T4

% cholesky: A = R'R
% Ax = b --> R'Rx = b
% y = Rx --> R'y = b

n = 4;
A = diag(5:8) - diag(ones(n-1, 1), 1) - diag(ones(n-1, 1), -1);
b = ones(n, 1);
R = chol(A);
R(2, 3);
y = b \ R';
y(2);

%% T5

% A a valori reali, quadrata, invertibile
% 1. F, ogni minore principale di A deve essere non singolare
% 2. F, deve essere < 1 il raggio spettrale della matrice B
% 3. V, per A SDP G-S converge sempre
% 4. F, se A tridiagonale O convergono O divergono sia J che G-S (se
%   convergono allora G-S converge più velocmente)

%% T6

n = 3;
A = [5 -1 0
    1 3 0
    0 0 6];
eig(A);

% 1. F, avendo due autovalori di modulo uguale QR non converge

x0_1 = ones(n, 1);
nmax = 2;
tol = 1e-6;
[lambda1, x1, iter1] = eigpower(A, tol, nmax, x0_1);

% 2. V, verificando applicando il metodo

x0_2 = [1 1 0]';
nmax = 100;
tol = 1e-6;
[lambda2, x2, iter2] = eigpower(A, tol, nmax, x0_2);

% 3. F, notando che x0_2 è multiplo dell'autovettore associato a lambda_2 =
% 4 e quindi ci converge subito, senza convergere all'autovalore massimo
% lambda_1 = 6

gershcircles(A)

% 4. V, verificata graficamente

%% T7

a = 1.5;
b = 3.5;
toll = 1e-4;
f = @(x) x .* sin(x);
Df = @(x) sin(x) + x.*cos(x);

[xvect, it] = bisez(a, b, toll, Df);

% verifica:
% f = @(x) -x * sin(x);
% fminbnd(f, 1.5, 3.5);
% x = a:.01:b;
% plot(x, f(x))

%% T8

% verificando l'annullamento delle derivate otteniamo moltplicità 3
f = @(x) (exp(x) - 1) * x^2;
df = @(x) exp(x)*x^2 + 2*exp(x)*x - 2*x;
x0 = 1.5;
toll = 1e-6;
mol = 3;

[xvect, it] = newton(x0, nmax, toll, f, df, mol);
xvect(2);
xvect(3);

%% T9

% come si fa???

%% T10

x(1) = 10;
iter = 10;

for k = 1:iter+1
    x(k+1) = 0.5 * (x(k) + 19 / x(k));
end

err = abs(sqrt(19) - x);

k = 0:iter;
semilogy(k, err(k+1), k, 1./k, k, 1./k.^2)

% quadratico, metodo grafico
% come posso garantire???

% huh????????
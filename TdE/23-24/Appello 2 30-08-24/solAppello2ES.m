clear
clc
close all

n = 100;
A = 6*diag(ones(n, 1)) - 4*diag(ones(n-1, 1), 1) - 4*diag(ones(n-1, 1), -1) ...
    + diag(ones(n-2, 1), 2) + diag(ones(n-2, 1), -2);
x_ex = ones(n, 1);
b = A * x_ex;

% A SDP

%% 1

[L, U, P] = lu(A);
y = fwsub(L, P * b);
x_hat = bksub(U, y);

% nflops solo per fattorizzazione
nflops_LU = 2/3 * n^3

spy(P)
isequal(P, eye(n));
% viene effettuata permutazione righe

% residuo
res = b - A * x_hat;
% residuo relativo (normalizzato)
res_norm = norm(res) / norm(b);
% Stima errore relativo
stima_err_rel = res_norm * cond(A)

% errore relativo
err_rel = norm(x_ex - x_hat) / norm(x_ex)

% essendo A SDP sarebbe convenuto fare cholesky
nflops_CHOL = 1/3 * n^3

%% 2

D = diag(diag(A));
B_j = eye(n) - D\A;
rho_j = max(abs(eig(B_j)))
% jacobi non converge

T = tril(A);
B_gs = eye(n) - T\A;
rho_gs = max(abs(eig(B_gs)))
if rho_gs < 1
    disp('GS converge')
end

x0 = zeros(n, 1);
toll = 1e-3;
nmax = 1e3;
[x_gs, k_gs] = gs(A,b,x0,toll,nmax);
k_gs

res = b - A * x_gs;
res_norm = norm(res) / norm(b)

err_rel = norm(x_ex - x_gs) / norm(x_ex)

% accuratezza molto scarsa

%% 3

P1 = eye(n);
P2 = 2*diag(ones(n,1)) - diag(ones(n-1,1),1) - diag(ones(n-1,1),-1);
P3 = 6*diag(ones(n, 1)) - 4*diag(ones(n-1, 1), 1) - 4*diag(ones(n-1, 1), -1);

eig(P3)
% P3 non SPD, non garantisce convergenza

K1 = max(eig(P1\A)) / min(eig(P1\A));
c1 = (sqrt(K1) - 1) / (sqrt(K1) + 1)
K2 = max(eig(P2\A)) / min(eig(P2\A));
c2 = (sqrt(K2) - 1) / (sqrt(K2) + 1)
% converge prima con P2

tol = 1e-6;
maxit = 100;
% x0 = zeros(n, 1) default
x = pcg(A,b,tol,maxit,P2);

err_norm_A = sqrt((x_ex - x)' * A * (x_ex - x))

%% 4

x0 = (1:n)';
y0 = x0 / norm(x0);
tol = 1e-3;
nmax = 1000;

[lambda,x,iter] = eigpower(A,tol,nmax,y0);
lambda

eigA = sort(eig(A), "descend");

% dalla teoria, per A reale simmetrica:
% |lambdaA(1) - lambda_k| <= C |(eigA(2) - /eigA(1)|^2k
% dato che se k = 0 vale:
% |lambdaA(1) - lambda_0| <= C

% la richiesta diventa dunque:
% |lambdaA(1) - lambda_k| / |lambdaA(1) - lambda_0| < 10^-3
% usando la maggiorazione precedente:

% |(eigA(2)/eigA(1)|^2k < 10^-3
% 2k * log(|(eigA(2)/eigA(1)|) < log(1e-3)
% dato che log(|(eigA(2)/eigA(1)|) < 0 cambio verso alla disequazione

N = 0.5 * log(1e-3) / log(abs(eigA(2)/eigA(1)));
ceil(N)

% dalla definizione di ordine di convergenza, il metodo delle potenze
% converge sempre con ordine 1
% p = 1

% per verifica
[lambdaV,x,iter] = eigpower_it(A,tol,nmax,y0);
[p,c] = stimap(lambdaV);

%% 5

x0 = zeros(n, 1);
mu = 1.7;
tol = 1e-3;

[x_vec, it] = es5_metodo_iterativo(A, b, x0, mu, tol);
it
x_vec(1, 2)
x_vec(1, 3)
x_vec(1, end)
clear
clc
close all

n = 100;
A = 21*diag(ones(n, 1)) - 20*diag(ones(n-1, 1), 1) - 4*diag(ones(n-1, 1), -1) ...
    + 2*diag(ones(n-2, 1), 2) + diag(ones(n-2, 1), -2);

x_ex = 2*ones(n, 1);
b = A * x_ex;

%% 1

P_J = diag(diag(A));
P_GS = tril(A);

B_J = eye(n) - P_J\A;
B_GS = eye(n) - P_GS\A;

rho_B_J = max(abs(eig(B_J)))
rho_B_GS = max(abs(eig(B_GS)))
% converge solo GS

%% 2

tol = 1e-3;
nmax = 1e4;
x0 = b;

[x,k] = gs(A,b,x0,tol,nmax);
% alternativamente chiamando richardson con alpha = 1 e P = P_GS
x(1)
k

res = b - A*x;
res_norm = norm(res) / norm(b)

%% 3

cond(A)
err_rel_est = cond(A) * res_norm

%% 4

T = tril(A);
I = eye(n);

% calcolo a mano B

% x_k+1 = B * x_k + g
% 1/omega * T * x_k+1 = - (A * + (1/omega - 2) * T) * x_k + b
% x_k+1 = - omega*T^-1 * (A * + (1/omega - 2) * T) * x_k + omega*T^-1 * b
% B = - omega*T^-1 * (A * + (1/omega - 2) * T)
% B = - omega * T^-1 * A - (1 - 2*omega) * I
% B = (2*omega - 1) * I - omega * T^-1 * A

% risolvo per via grafica per trovare il valore di omega per cui troviamo
% il raggio spettrale minimo

rho_B_vec = [];
omega_vec = 0:0.01:2;
for omega = omega_vec
    B = (2*omega - 1) * I - omega * (T\A);
    rho_B = max(abs(eig(B)));
    rho_B_vec = [rho_B_vec, rho_B];
end

plot(omega_vec, rho_B_vec, omega_vec, ones(size(omega_vec)))

% zoom sull'intersezione con 1
clear rho_B_vec
rho_B_vec = [];
omega_vec = 1.25:0.0001:1.35;
for omega = omega_vec
    B = (2*omega - 1) * I - omega * (T\A);
    rho_B = max(abs(eig(B)));
    rho_B_vec = [rho_B_vec, rho_B];
end
figure
subplot(1,2,1)
plot(omega_vec, rho_B_vec, omega_vec, ones(size(omega_vec)))

% zoom sul minimo
clear rho_B_vec
rho_B_vec = [];
omega_vec = 0.75:0.0001:0.85;
for omega = omega_vec
    B = (2*omega - 1) * I - omega * (T\A);
    rho_B = max(abs(eig(B)));
    rho_B_vec = [rho_B_vec, rho_B];
end
subplot(1,2,2)
plot(omega_vec, rho_B_vec, omega_vec, ones(size(omega_vec)))

% dallo zoom
[M, I] = min(rho_B_vec);
M
omega_vec(I)

%% 5

% la matrice non è SDP

AtA = A' * A;

tol = 1e-6;
nmax = 1e3;
x0 = ones(n, 1);

% calcolo autovalore minimo e massimo
[lambda,x,iter] = eigpower(AtA,tol,nmax,x0);
[mu,x,iter] = invpower(AtA,tol,nmax,x0);

% calcolo il numero di condizionamento in norma 2
K2 = sqrt(lambda / mu)

% verifico
cond(A)

%% 6

% uso la fattorizzazione QR per trovare la soluzione ai minimi quadrati
C = A(:, 1:end-2);
[Q, R] = qr(C);
z = R\(Q' * b);
z(1)

%% 7

% copio la funzione vettoriale (sistema di equazioni)
F = @(x) A*x + exp(-x/10) - b - ones(n, 1);

x0 = b;
B0 = A;
tol = 1e-4;
nmax = 100;

% chiamo il metodo che ho scritto
[xvect, it] = bfgs_zero(F, B0, x0, tol, nmax);
xvect(1, 2:4)
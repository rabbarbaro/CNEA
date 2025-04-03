clear
clc
close all

n = 100;
A = 20*diag(ones(n, 1)) - 11*diag(ones(n-1, 1), 1) - 11*diag(ones(n-1, 1), -1) ...
    + diag(ones(n-2, 1), 2) + diag(ones(n-2, 1), -2);

%% 1

P_J = diag(diag(A));
P_GS = tril(A);

B_J = eye(n) - P_J\A;
B_GS = eye(n) - P_GS\A;

rho_B_J = max(abs(eig(B_J)))
rho_B_GS = max(abs(eig(B_GS)))
% J non converge
% GS converge più velocemente di J

% fatt_abb_err = rho_B_GS^k

%% 2

b = 5*ones(n, 1);
tol = 1e-2;
nmax = 1e4;
x0 = b;

[x,k] = gs(A,b,x0,tol,nmax);
x(1)
k

res = b - A*x;
res_norm = norm(res) / norm(b)

%% 3

cond(A)
err_rel_est = cond(A) * res_norm

%% 4

% P = 1/omega * tril(A);
% B = eye(n) - P\A = I - omega * T^-1 * A
% g = P\b = omega * T^-1 * b;

rho_B_vec = [];
for omega = [1.45, 1.55, 1.65, 1.75, 1.85]
    P = 1/omega * tril(A);
    B = eye(n) - P\A;
    rho_B = max(abs(eig(B)));
    rho_B_vec = [rho_B_vec, rho_B];
end

rho_B_vec(3)

%% 5

P1 = eye(n);
P2 = 2*diag(ones(n, 1)) - diag(ones(n-1, 1), 1) - diag(ones(n-1, 1), -1);

K_P1A = max(eig(P1\A)) / min(eig(P1\A));
K_P2A = max(eig(P2\A)) / min(eig(P2\A))
% converge più rapidamente con P2

k = 10;
d_P2 = (K_P2A - 1) / (K_P2A + 1);
fatt_abb_P2_grad = d_P2^k

%% 6

c_P2 = (sqrt(K_P2A) - 1) / (sqrt(K_P2A) + 1);
fatt_abb_P2_conj_grad = (2 * c_P2^k) / (1 + c_P2^(2*k))

%% 7

% fatto a esercitazione
% K = sdpcond_mod(A, 1e-12, 1)
% K = sdpcond_mod(A, 1e-12, 2)
% K = sdpcond_mod(A, 1e-12, 100)

K = K_SDP(A, 1e-12, 1, ones(n, 1))
K = K_SDP(A, 1e-12, 2, ones(n, 1))
K = K_SDP(A, 1e-12, 100, ones(n, 1))

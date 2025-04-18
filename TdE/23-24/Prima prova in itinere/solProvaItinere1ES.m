clear
clc
close all

n = 10;
A = 8*diag(ones(n, 1)) - 5*diag(ones(n-1, 1), 1) - 5*diag(ones(n-1, 1), -1) ...
    + diag(ones(n-2, 1), 2) + diag(ones(n-2, 1), -2);
x_ex = ones(n, 1);
b = A * x_ex;

%% 1

% A è SDP, quindi sicuramente GS converge per ogni x0

D_inv = diag(1 ./ diag(A));
B_j = eye(n) - D_inv * A;       % oppure B_j = eye(n) - diag(diag(A))\A;
rho_B_j = max(abs(eig(B_j)));
% è > 1, non converge

P_gs = tril(A);
B_gs = eye(n) - P_gs \ A;
rho_B_gs = max(abs(eig(B_gs)));
% verificato che converge

x0 = zeros(n, 1);
toll = 1e-6;
nmax = 1000;
[x, k] = gs(A, b, x0, toll, nmax);

res_norm = norm(b - A*x) / norm(b);

err_rel_est = cond(A) * res_norm;
err_rel = norm(x_ex - x) / norm(x_ex);

%% 2

P1 = eye(n);
P2 = 9*diag(ones(n, 1)) - 4*diag(ones(n-1, 1), 1) - 4*diag(ones(n-1, 1), -1);
P3 = 2*diag(ones(n, 1)) - diag(ones(n-1, 1), 1) - diag(ones(n-1, 1), -1);

max(eig(P1\A))/min(eig(P1\A))
max(eig(P2\A))/min(eig(P2\A))
max(eig(P3\A))/min(eig(P3\A))
% P3

x0 = zeros(n, 1);
tol = 1e-6;

[x, k] = richardson(A, b, P3, x0, tol, nmax);

err = x_ex - x;
err_norm_A = sqrt(err' * A * err)

%% 3

% it = n = 10

maxit = 10;
x0 = zeros(n, 1);
tol = 1e-6;

x = pcg(A,b,tol,maxit,eye(n),eye(n),x0)

err_pcg = x_ex - x;
norm(err_pcg)
err_norm_A = sqrt(err_pcg' * A * err_pcg)
% ??????????

%% 4

D = qrbasic(A,tol,3);
s = D(9);
x0 = ones(n, 1);
x0(1:2:end) = -1;
[lambda,x,iter] = invpowershift(A,s,1e-3,nmax,x0)

x_norm = x / norm(x)

%% 5

x0 = ones(n, 1);
xV = [];
tol = 1e-3;
gamma = 0.1;
phi = @(y) sin(1/n * y' * y) + 0.5 * y' * A * y;
grad_phi = @(y) 2/n * cos(1/n * y' * y) * y + A * y;

k = 0;
while norm(grad_phi(x0)) > tol
    x = x0 - gamma * grad_phi(x0);
    xV = [xV x];
    gamma = abs((x - x0)' * (grad_phi(x) - grad_phi(x0))) / norm(grad_phi(x) - grad_phi(x0))^2;
    x0 = x;
    k = k+1;
end
xV
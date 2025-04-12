clear
clc
close all

n = 100;
A = 151/5 * diag(ones(n,1)) -16 * diag(ones(n-1,1),1) -16 * diag(ones(n-1,1),-1)...
    + diag(ones(n-2,1),2) + diag(ones(n-2,1),-2);
vec = (1:n)';
b = (n + vec)/n;

% A SPD

%% 1

% applico gs
x0 = zeros(n, 1);
toll = 1e-6;
nmax = 1000;
[x_gs, Nit] = gs(A,b,x0,toll,nmax);
% Nit = 978

% Residuo
res = b - A * x_gs;
% Residuo relativo/normalizzato
res_norm = norm(res) / norm(b);
% res_norm = 9.9605e-07

% Stima (maggiorazione) errore relativo (quando deltaA ~ 0, wilkinson-higham)
stima_err_rel = res_norm * cond(A)
% stima_err_rel = 3.0205e-04

% risolvo per via grafica, converge se rho_B < 1
% zoommo su [0.85, 0.9]
rho_B_v = [];
beta_v = 0.85:0.0001:0.9;
for beta = beta_v
    D = diag(diag(A));
    P = beta * D + (1 - beta) * tril(A);
    % Metodo iterativo
    B = eye(n) - P \ A;
    rho_B = max(abs(eig(B)));
    rho_B_v = [rho_B_v rho_B]
end

plot(beta_v, rho_B_v, beta_v, ones(length(beta_v)))
% circa 0.8946

%% 2

% gradiente coniugato precondizionato

% definisco i precondizionatori
P1 = 151/5 * eye(n);
P2 = 2*diag(ones(n,1)) - diag(ones(n-1,1),1) - diag(ones(n-1,1),-1);
P3 = 3*diag(ones(n,1)) - diag(ones(n-1,1),1) - diag(ones(n-1,1),-1);
% sono sdp

% calcolo i numeri di condizionamento spettrali e i parametri c
% quella con il parametro c più piccolo converge prima
K_invP1A = max(eig(P1\A)) / min(eig(P1\A));
c1 = (sqrt(K_invP1A) - 1) / (sqrt(K_invP1A) + 1)
% c1 = 0.8914
K_invP2A = max(eig(P2\A)) / min(eig(P2\A));
c2 = (sqrt(K_invP2A) - 1) / (sqrt(K_invP2A) + 1)
% c2 = 0.6093
K_invP3A = max(eig(P3\A)) / min(eig(P3\A));
c3 = (sqrt(K_invP3A) - 1) / (sqrt(K_invP3A) + 1)
% c3 = 0.7725

% converge prima con P2

% vale ||x-xk||_A <= (2*c^k) / (1 + c^(2k)) * ||x-x0||_A;
% dato che per k grandi vale circa
% f_abb_k = 2 * c^k;
% % Numero minimo di iterazioni data tolleranza

tol = 1e-6;
k_min = ceil(log(tol / 2) / log(c2))

x0 = zeros(n, 1);
[x, ~, rel_res] = pcg(A,b,tol,1000,P2,eye(n),x0);
% pcg converged at iteration 9 to a solution with relative residual 7.2e-07
rel_res

%% 3

% applico il metodo
x0 = ones(n, 1);
nmax = 1000;
tol = 1e-6;
[lambda,x,iter] = invpower(A,tol,nmax,x0);
% Nit = 15

% autovalore zero
y = x0 / norm(x0);
lambda0 = y' * A * y
% autovalori 1,2
[lambda1,~,iter] = invpower(A,tol,1,x0)
[lambda2,~,iter] = invpower(A,tol,2,x0)
% lambda0 =  0.4800
% lambda1 = 0.2199
% lambda2 = 0.2144
% lambda = 0.2116 dopo Nit = 15

%% 4

n = 1000;
x0 = ones(n,1);
tol = 1e-2;
nmax = 1000;

[lambda_vec, x_seq, Nit] = grad_sfera(A, tol, nmax, x0);

lambda0 = lambda_vec(1)
lambda1 = lambda_vec(2)
lambdaN = lambda_vec(end)

function [lambda_vec, xvec, N_iter] = grad_sfera(A, tol, nmax, x0)

end

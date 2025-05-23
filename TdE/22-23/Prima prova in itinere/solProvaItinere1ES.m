clear
clc
close all

n = 225;
A = full( gallery('poisson', 15 ));
b = 2*ones(n, 1);

% A SDP

%% 1

% verifico se è simmetrica
if isequal(A, A')
    disp('A simmetrica');
    nA = size(A, 1);

    % calcolo il determinante delle sottomatrici (criterio di Sylvester)
    for ii = 1:nA
        if det(A(1:ii, 1:ii)) > 0
            % se siamo arrivati all'ultima iterazione allora le condizioni
            % su tutti i determinanti sono verificate
            if ii == nA
                disp('A DP');
            end
        else
            error('A non DP');
        end
    end
else
    disp('A non simmetrica');
end

% cholesky
est_flops = 1/3 * n^3 + 2 * n^2

%% 2

R = chol(A);        % R triangolare superiore
y = R' \ b;         % R'y = b
x_cap = R \ y;      % Rx = y

res = b - A*x_cap;
res_rel = norm(res) / norm(b)
err_rel_est = cond(A) * res_rel

%% 3

D_inv = diag(1 ./ diag(A));
B_j = eye(n) - D_inv * A;       % oppure B_j = eye(n) - diag(diag(A))\A;
rho_B_j = max(abs(eig(B_j)))

k = 100;
fatt_abb = rho_B_j^k

%% 4

x0 = b;
toll = 1e-12;
nmax = 1;
[xk,k] = richardson(A, b, eye(n), x0, toll, nmax);

phi = @(y) 1/2 * y' * A * y - y' * b;
phi(x0)
phi(xk)

%% 5

P1 = 4*diag(ones(n, 1)) - diag(ones(n-1, 1), 1) - diag(ones(n-1, 1), -1);
R2 = ichol(sparse(A));
P2 = R2' * R2;

K1 = max(eig(P1\A)) / min(eig(P1\A))
K2 = max(eig(P2\A)) / min(eig(P2\A))

% P2

k = 20;
d2 = (K2 - 1) / (K2 + 1);
fatt_abb = d2^k

%% 6

x0 = ones(n, 1);
nmax = 1000;
tol = 1e-6;
[lambda, x, iter] = eigpower_it(A, tol, nmax, x0);

[p,c] = stimap(lambda);

%% 7

tol = 1e-6;
D = qrbasic_shift(A,tol,1);
D(2)
D = qrbasic_shift(A,tol,2);
D(2)
D = qrbasic_shift(A,tol,20);
D(2)

D = qrbasic(A, tol, 1000);

eig(A);